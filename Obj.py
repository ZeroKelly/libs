from pyecharts.charts import Map, Scatter, Pie
from pyecharts import options as opts
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import pandas as pd
import json
import math
import re
import threading
import time
import hail as hl

from genodig import show
from genodig.common import uid
from genodig import geno
from genodig import pheno
from genodig.fields import fields

composition_data = pd.read_csv('s3://genodig-prod/workspaces/users/ancestry/basic_data/民族成分key.csv')
composition_key = composition_data.val.tolist()
composition_key.remove('唾液盒编号')
composition_key.remove('东亚人')

del composition_data

def valid_mean(series):
    high = series.mean()+1.5*series.std()
    low = series.mean()-1.5*series.std()
    return series[(series<=high) & (series>=low)].mean()


def map_chip_version(L):
    if L is None:
        return None
    elif len(L) < 1:
        return None
    else:
        chip_version_map = {'v2.1':'affy-v2.1', 'v2.2':'affy-v2.2', 'v2':'affy-2', 'v1':'affy-1', 'i1':'ilmn-v1.0', 'i1.1':'ilmn-v1.1'}
        return [chip_version_map[i] for i in L]
    
def add_weighted_hap_ratio_info(df, weighted_hap_ratio_data):
    weighted_hap_ratio_data = weighted_hap_ratio_data[['yhaplogroup', 'ratio']]
    weighted_hap_ratio_data.columns = ['yhaplogroup', 'total_ratio']
    df = pd.merge(df, weighted_hap_ratio_data, on='yhaplogroup', how='left')
    df['是该单倍群比例均值的倍数'] = df['ratio']/df['total_ratio']
    df['是该单倍群比例均值的倍数'][df['是该单倍群比例均值的倍数'] == pd.np.inf] = 0
    df.sort_values('是该单倍群比例均值的倍数', ascending=False, inplace=True)
    return df


def fix_data(original_data, fix_data, merge_key, fix_val_key):
    bad_data = original_data[original_data[merge_key].isin(fix_data[merge_key])]
    good_data = original_data[original_data[merge_key].isin(fix_data[merge_key]) == False]
    new_key = bad_data.keys().drop(fix_val_key)
    bad_data = bad_data[new_key]
    fix = pd.merge(bad_data, fix_data, on=merge_key, how='left')
    fix = fix[good_data.keys()]
    new_data = fix.append(good_data)
    return new_data

class Data:
    def __init__(self, data, tree):
        self.full_data = data
        self.data = data
        self.tree = tree
        self.map = None
        self.select = None
    
    def use_select_data(self):
        self.tmp_data = self.data.copy(deep=True)
        self.data = self.select
        
    def unuse_select_data(self):
        self.data = self.tmp_data
       
    def chip_version_distribution(self, hapname):
        freq = self.full_data[self.full_data['yhaplogroup'] == hapname].version.value_counts()
        return freq
    
    def add_isogg_map(self, isogg_map):
        self.map = isogg_map
    
    def generate_isoggname_map(self):
        tmp = self.full_data[['yhaplogroup', 'y23andmename']]
        self.map = tmp.drop_duplicates('yhaplogroup')
        self.map.columns = ['yhaplogroup', 'isogg']
        
    def add_isogg_name(self, df):
        if self.map is None:
            self.generate_isoggname_map()
        df = pd.merge(df, self.map, on='yhaplogroup', how='left')
        new_order = ['yhaplogroup', 'isogg']+ df.keys().drop(['yhaplogroup', 'isogg']).tolist()
        return df[new_order]
    
    def hap_ratio_by_surname(self, surname):
        surname_num = self.select_by_surname(surname)
        freq = surname_num['yhaplogroup'].value_counts()
        res = pd.DataFrame()
        res['yhaplogroup'] = freq.keys().tolist()
        res['num'] = freq.values
        res['ratio'] = res['num']/surname_num.shape[0]
        res.sort_values('ratio', ascending=False, inplace=True)
        res = self.add_isogg_name(res)
        return res
    
    def hap_ratio_by_nation(self, nation):
        surname_num = self.select_by_nation(nation)
        freq = surname_num['yhaplogroup'].value_counts()
        res = pd.DataFrame()
        res['yhaplogroup'] = freq.keys().tolist()
        res['num'] = freq.values
        res['ratio'] = res['num']/surname_num.shape[0]
        res.sort_values('ratio', ascending=False, inplace=True)
        res = self.add_isogg_name(res)
        return res
    
    def initialize_selector(self):
        self.select = self.data
        print('初始化选择器，包含 %s 条数据'%(self.select.shape[0]))
    
    def selector(self, **kwargs):
        tmp = self.select.copy(deep=True)
        for i in kwargs:
            if type(kwargs[i]) == list:
                tmp = tmp[tmp[i].isin(kwargs[i])]
            else:
                tmp = tmp[tmp[i] == kwargs[i]]
        self.select = tmp
        print('已选中%s条数据'%(self.select.shape[0]))

    def composition_selector(self, **kwargs):
        tmp = self.select.copy(deep=True)
        for i in kwargs:
            tmp = tmp[tmp[i] >= kwargs[i]]
        self.select = tmp
        print('已选中%s条数据'%(self.select.shape[0]))
    
    def select_version(self, versions=None):
        if versions is None:
            self.data=self.full_data
            return
        elif len(versions) < 1:
            self.data=self.full_data
        else:
            self.data = self.full_data[self.full_data['version'].isin(versions)]
    
    def select_by_surname(self, surname):
        tmp = self.data[self.data['surname'] == surname]
        return tmp

    def select_by_hapname(self, hapname):
        tmp = self.data[self.data['yhaplogroup'] == hapname]
        return tmp
    
    def select_by_prov(self, prov):
        tmp = self.data[self.data['province'] == prov]
        return tmp
    
    def select_by_city(self, city):
        tmp = self.data[self.data['city'] == city]
        return tmp
    
    def select_by_nation(self, nation):
        tmp = self.data[self.data['nation'] == nation]
        return tmp
    
    def select_by_area(self, area):
        tmp = self.data[self.data['area'] == area]
        return tmp
    
    def select_by_uid(self, uid):
        return self.data[self.data['uid'].isin(uid)]
    
    def select_by_barcode(self, barcode):
        return self.data[self.data['barcode'].isin(barcode)]

    def select_by_hapnames(self, hapnames):
        return self.data[self.data['yhaplogroup'].isin(hapnames)]

    def select_haplogroup(self, hapname):
        '''select haplogroup and all subhaps'''
        return self.select_by_hapnames(self.tree.get_branches(hapname))

    def select_haplogroup_star(self, hapname):
        '''select only star hap'''
        return self.select_by_hapname(hapname)
    
    def select_branch(self, hapname, include_sub_hap=False):
        if include_sub_hap==False:
            return self.select_by_hapname(hapname)
        else:
            haps = self.tree.get_branches(hapname)
            return self.select_by_hapnames(haps)

    def cal_ratio_table(self, sub_df, total_df, key='province', debug=False):
        df1 = pd.DataFrame()
        freq = sub_df[key].value_counts()
        df1['name'] = freq.keys().tolist()
        df1['value'] = freq.tolist()
        df2 = pd.DataFrame()
        freq = total_df[key].value_counts()
        df2['name'] = freq.keys().tolist()
        df2['value'] = freq.tolist()
        result = pd.merge(df1, df2, on='name', how='inner')
        result['ratio'] = result.iloc[:, 1] * 1.0 / result.iloc[:, 2]
        if debug:
            return result
        result = result[['name', 'ratio']]
#         result.loc[result.shape[0] + 1] = ['mean', sub_df.shape[0] * 100.0 / total_df.shape[0]]
        result = result.sort_values('ratio', ascending=False)
        return result

    def hap_distribution(self, hapname, key='province'):
        sub_haps = self.tree.get_branches(hapname)
        tmp = self.select_by_hapnames(sub_haps)
        ratio = self.cal_ratio_table(tmp, self.data, key)
        ratio.columns = [key, hapname]
        return ratio

    def hap_distribution_star(self, hapname, key='province'):
        tmp = self.select_by_hapname(hapname)
        ratio = self.cal_ratio_table(tmp, self.data, key)
        ratio.columns = [key, hapname + '*']
        return ratio
    
#     def city_hap_ratio(self, hapname, include_sub_hap):
#         if include_sub_hap:
#             res = self.hap_distribution(hapname, 'city')
#         else:
#             res = self.hap_distribution_star(hapname, 'city')
#         res.columns = ['city', 'ratio']
#         return res
    
    def cal_province_hap_ratio(self, prov, hap, include_sub_hap, city_population_info):
        city_comp = self.data[self.data['province'] == prov].city.value_counts()

        cities = city_population_info[city_population_info['province'] == prov].city.unique().tolist()

        sub_df = self.select_branch(hap, include_sub_hap)
        hap_city_comp = sub_df[sub_df['province'] == prov].city.value_counts()
        
        pop = city_population_info[city_population_info.city.isin(cities)][['province','city', 'population']]

        sample_num = pd.DataFrame()
        hap_sample_num = pd.DataFrame()
        hap_sample_num['city'] = hap_city_comp.keys().tolist()
        hap_sample_num['hap_sample_num'] = hap_city_comp.values.tolist()
        sample_num['city'] = city_comp.keys().tolist()
        sample_num['sample_num'] = city_comp.values.tolist()

        res = pd.merge(hap_sample_num, sample_num, on='city', how='right')
        res = pd.merge(res, pop, on='city', how='right')

        res['hap_ratio'] = res['hap_sample_num']/res['sample_num']

        res['hap_pop'] = res['population']*res['hap_ratio']
        
        return res    
    
    
# #     算全国的时候也应该加权算？
#     def hap_ratio_eval(self, hap, include_sub_hap, province_population_info, debug=False):
#         tmp = self.distribution_over_provinces(hap, include_sub_hap)

#         province_population_info = province_population_info[['province', 'population']]

#         res = pd.merge(tmp, province_population_info, on='province', how='outer')

#         res['hap_pop'] = res['population']*res['ratio']

#         hap_ratio_eval = res.hap_pop.sum()/res.population.sum()
        
#         if debug:
#             print('%s / %s = %s'%(res.hap_pop.sum(), res.population.sum(), hap_ratio_eval))
#             return res

#         return hap_ratio_eval
    
    # 市级加权->省->全国
    def hap_ratio_eval(self, hap, include_sub_hap, city_population_info):
        tmp = self.distribution_over_city(hap, include_sub_hap, city_population_info)
        eval_ratio = tmp['hap_pop'].sum()/tmp['population'].sum()
        return eval_ratio
        
        
    def distribution_over_city(self, hap, include_sub_hap, city_population_info):
        provs = city_population_info.province.unique()
        result = pd.DataFrame()
        for i in provs:
            result = result.append(self.cal_province_hap_ratio(i, hap, include_sub_hap, city_population_info))
        result = result[pd.isna(result['sample_num'])==False]
        return result
    
    def distribution_in_province(self, prov, hap, include_sub_hap, city_population_info):
        provs = self.distribution_over_city(hap, include_sub_hap, city_population_info)
        provs = provs[provs['province'] == prov]
        return provs
    
    def distribution_over_provinces_kaijun_method(self, hap, include_sub_hap, city_population_info):
        provs = city_population_info.province.unique()
        result = self.distribution_over_city(hap, include_sub_hap, city_population_info)
        res = []
        for i in provs:
            tmp = result[result['province'] == i]
            val = tmp['hap_pop'].sum()/tmp['population'].sum()
            res.append([i, val])
        res = pd.DataFrame(res, columns=['province', 'ratio'])
        res.sort_values('ratio', ascending=False, inplace=True)
        return res
    
    def distribution_over_provinces(self, hapname, include_sub_hap=True):
        tmp = self.select_branch(hapname, include_sub_hap)
        tb = self.cal_ratio_table(tmp, self.data, key='province')
        tb.rename(columns={'name':'province'}, inplace=True)
        return tb

    def decompose_direct_branch(self, hapname, key='province'):
        total = self.hap_distribution(hapname)
        star = self.hap_distribution_star(hapname, key)
        total = pd.merge(total, star, on=key, how='outer')
        for i in self.tree.get_direct_branches(hapname):
            sub_hap = self.hap_distribution(i)
            total = pd.merge(total, sub_hap, on=key, how='outer')
        #         total = self.formalization_tb(total, key)
        return total

    def exclude_hap(self, hap_list, exclude_hap):
        exclude_haps = set(self.tree.get_branches(exclude_hap))
        return list(set(hap_list).difference(exclude_haps))

    def exclude_haps(self, hap_list, exclude_hap_list):
        new_hap_list = set(hap_list[:])
        for i in exclude_hap_list:
            new_hap_list = new_hap_list.difference(self.tree.get_branches(i))
        return list(new_hap_list)

    def customize_decompose(self, hap_list, key='province'):
        total = self.hap_distribution(hap_list[0], key)
        for i in hap_list[1:]:
            if '*' in i:
                sub_hap = self.hap_distribution_star(i[:-1])
            else:
                sub_hap = self.hap_distribution(i)
            total = pd.merge(total, sub_hap, on=key, how='outer')
        total['remain'] = total[hap_list[0]]
        for i in hap_list[1:]:
            total['remain'] = total['remain'] - total[i]
        total.sort_values('remain', ascending=False, inplace=True)
        return total

    def decompse_haplogroup(self, hapname, key='province'):
        total = self.hap_distribution(hapname)
        for i in self.tree.get_branches(hapname):
            sub_hap = self.hap_distribution_star(i)
            total = pd.merge(total, sub_hap, on=key, how='outer')
        #         total = self.formalization_tb(total, key)
        return total

    def append_haps(self, haps):
        total = haps[0]
        for i in haps:
            total = total.append(i)
        return total

    def distribution_analyze(self, hapname, key='province'):
        total = self.hap_distribution(hapname)
        star = self.hap_distribution_star(hapname, key)

    def cal_ratio_to_top(self, tb, key='province'):
#         total = tb[tb[key] == 'mean'].iloc[0, 1]
#         tb.loc[tb.shape[0] + 1] = ['ratio_to_top'] + (tb[tb['province'] == 'mean'].iloc[:, 1:] * 1.0 / total).iloc[
#             0].tolist()
        return tb

    def cal_bc(self, tb):
        haps = [tb.keys()[1]] * 2 + tb.keys()[3:].tolist()
        times = []
        for i in haps:
            times.append(self.tree.get_bc_time(i))
        tb.loc[tb.shape[0] + 1] = ['time(BC)'] + times
        return tb
    
    def cal_diversity(self, tb):
        tb['diversity'] = ((tb.iloc[:, 2:]>0)*1).sum(axis=1)
        tb.sort_values(['diversity', tb.keys()[1]], ascending=False, inplace=True)
        return tb

    def formalization_tb(self, tb, key='province'):
#         tb = self.cal_ratio_to_top(tb, key)
        tb = self.cal_bc(tb)
        tb.sort_values(tb.keys()[1], ascending=False, inplace=True)
        return tb
    
    def draw_distribution_map(self, hapname, include_sub_hap=False):
        tb = self.distribution_over_provinces(hapname, include_sub_hap)
        max_val, min_val = tb.iloc[:, 1].max(), tb.iloc[:, 1].min()
        map_data = []
        for i in range(tb.shape[0]):
            tmp = tb.iloc[i, :]
            map_data.append([tmp.iloc[0], tmp.iloc[1]])
        DIC = {True:'（含下游）', False:'（不含下游）'}
        m = Map()
        m.add('', map_data, "china", is_map_symbol_show=False, is_selected=False)
        m.set_global_opts(
            title_opts=opts.TitleOpts(title=hapname+' '+DIC[include_sub_hap]+'分布频率'),
            visualmap_opts=opts.VisualMapOpts(
                is_show=False,
                orient='horizontal',\
                max_=max_val, min_=min_val, \
                range_color=['#FFFFFF','#C4DCF9','#99CAFF','#66B0FF','#3395FF'],\
                range_text = ['占比较大', '占比较小'],
#                 split_number=5,\
#                 is_piecewise=True
            ),
            legend_opts=opts.LegendOpts(is_show=False)
        )
        m.set_series_opts(
            label_opts=opts.LabelOpts(is_show=False)
        )
#         mean_val = ['mean']+tb.mean().tolist()
        sum_val = ['sum']+[self.hap_ratio_eval(hapname, include_sub_hap, city_population_info)]
#         tb.loc[tb.shape[0]] = mean_val
        tb.sort_values('ratio', ascending=False, inplace=True)
        tb.loc[tb.shape[0]] = sum_val
        tb.rename(columns={'ratio':hapname+DIC[include_sub_hap]+'人口在各省占比', 'ratio1':'该省'+hapname+DIC[include_sub_hap]+'人口在全国占比'}, inplace=True)
        return m, tb
    
    def draw_distribution_map_weighted_method(self, hapname, include_sub_hap, city_population_info):
        tb = self.distribution_over_provinces_kaijun_method(hapname, include_sub_hap, city_population_info)
        max_val, min_val = tb.iloc[:, 1].max(), tb.iloc[:, 1].min()
        map_data = []
        for i in range(tb.shape[0]):
            tmp = tb.iloc[i, :]
            map_data.append([tmp.iloc[0], tmp.iloc[1]])
        DIC = {True:'（含下游）', False:'（不含下游）'}
        m = Map()
        m.add('', map_data, "china", is_map_symbol_show=False, is_selected=False)
        m.set_global_opts(
            title_opts=opts.TitleOpts(title=hapname+' '+DIC[include_sub_hap]+'分布频率'),
            visualmap_opts=opts.VisualMapOpts(
                is_show=False,
                orient='horizontal',\
                max_=max_val, min_=min_val, \
                range_color=['#FFFFFF','#C4DCF9','#99CAFF','#66B0FF','#3395FF'],\
                range_text = ['占比较大', '占比较小'],
#                 split_number=5,\
#                 is_piecewise=True
            ),
            legend_opts=opts.LegendOpts(is_show=False)
        )
        m.set_series_opts(
            label_opts=opts.LabelOpts(is_show=False)
        )
#         mean_val = ['mean']+tb.mean().tolist()
        sum_val = ['sum']+[self.hap_ratio_eval(hapname, include_sub_hap, city_population_info)]
#         tb.loc[tb.shape[0]] = mean_val
        tb.sort_values('ratio', ascending=False, inplace=True)
        tb.loc[tb.shape[0]] = sum_val
        tb.rename(columns={'ratio':hapname+DIC[include_sub_hap]+'人口在各省占比', 'ratio1':'该省'+hapname+DIC[include_sub_hap]+'人口在全国占比'}, inplace=True)
        return m, tb
        
    def draw_distribution_in_prov_map_weighted_method(self, province, hapname, include_sub_hap, city_population_info):
        tb = self.distribution_in_province(province, hapname, include_sub_hap, city_population_info)[['city', 'hap_ratio']]
        tb.columns = ['province', 'ratio']
        max_val, min_val = tb.iloc[:, 1].max(), tb.iloc[:, 1].min()
        map_data = []
        for i in range(tb.shape[0]):
            tmp = tb.iloc[i, :]
            map_data.append([tmp.iloc[0], tmp.iloc[1]])
        DIC = {True:'（含下游）', False:'（不含下游）'}
        m = Map()
        m.add('', map_data, province, is_map_symbol_show=False, is_selected=False)
        m.set_global_opts(
            title_opts=opts.TitleOpts(title=hapname+' '+DIC[include_sub_hap]+'分布频率'),
            visualmap_opts=opts.VisualMapOpts(
                is_show=False,
                orient='horizontal',\
                max_=max_val, min_=min_val, \
                range_color=['#FFFFFF','#C4DCF9','#99CAFF','#66B0FF','#3395FF'],\
                range_text = ['占比较大', '占比较小'],
#                 split_number=5,\
#                 is_piecewise=True
            ),
            legend_opts=opts.LegendOpts(is_show=False)
        )
        m.set_series_opts(
            label_opts=opts.LabelOpts(is_show=False)
        )
#         mean_val = ['mean']+tb.mean().tolist()
#         sum_val = ['sum']+[self.hap_ratio_eval(hapname, include_sub_hap, city_population_info)]
#         tb.loc[tb.shape[0]] = mean_val
        tb.sort_values('ratio', ascending=False, inplace=True)
#         tb.loc[tb.shape[0]] = sum_val
        tb.rename(columns={'province':'city', 'ratio':hapname+DIC[include_sub_hap]+'人口在%s省占比'%(province), 'ratio1':'该省'+hapname+DIC[include_sub_hap]+'人口在全国占比'}, inplace=True)
        return m, tb
   

    def draw_selected_data_in_province(self, province):
        map_data = []
        freq = self.select.city.value_counts()
        for i in freq.keys():
            map_data.append([i, int(freq[i])])
        m = Map()
        m.add('', data_pair=map_data, maptype=province, is_map_symbol_show=False, is_selected=False)
        m.set_global_opts(
            visualmap_opts=opts.VisualMapOpts(
                is_show=False,
                orient='horizontal',\
                max_=max_val, min_=min_val, \
                range_color=['#FFFFFF','#C4DCF9','#99CAFF','#66B0FF','#3395FF'],\
                range_text = ['占比较大', '占比较小'],
#                 split_number=5,\
#                 is_piecewise=True
            ),
            legend_opts=opts.LegendOpts(is_show=False)
        )
        m.set_series_opts(
            label_opts=opts.LabelOpts(is_show=False)
        )
        return m

    def cal_distribution_center(self, haplogroup, include_sub_hap, city_population, coord, keep_ratio=0.1):
        tmp = self.distribution_over_city(haplogroup, include_sub_hap, city_population)
        tmp.sort_values('hap_ratio', ascending=False, inplace=True)
        tmp['coord'] = tmp.city.map(lambda x: coord.get(x, None))
        tmp.dropna(how='any', inplace=True)
        tmp['coordx'] = tmp.coord.map(lambda x: x[1])
        tmp['coordy'] = tmp.coord.map(lambda x: x[0])
        keep_data = tmp.iloc[:int(tmp.shape[0]*keep_ratio), :]
        x, y = valid_mean(keep_data.coordx), valid_mean(keep_data.coordy)
        return (x,y)
    
    def draw_migration_path(self, haplogroup, limit_year, city_population, direct_branch_count_data):
        tmp = self.tree.select_path(haplogroup, limit_year, direct_branch_count_data)
        for i in range(tmp.shape[0]):
            row = tmp.iloc[i, :]
            hap = row.yhaplogroup
            ybp = row.ybp
            branch_num = row.branch_num
            m,tb = self.draw_distribution_map_weighted_method(hap, False, city_population)
            if tb.iloc[:, 1].sum() == 0:
                continue
            else:
                print('%s ybp  %s direct_branches'%(ybp, branch_num))
                display(m.render_notebook())
    
    def draw_composition_pie_chart(self, data):
        res = data[composition_key].sum()
        res = res/data.shape[0]
        res = res[res>0]
        plt_data = [list(i) for i in zip(res.keys(), res.values)]
        m = Pie()
        m.add('', plt_data)
        m.set_global_opts(title_opts=opts.TitleOpts(title="民族成分"), legend_opts = opts.LegendOpts(is_show=False))
        return (m, res)
    
    def draw_pie_chart(self, data):
        plt_data = [list(i) for i in zip(data.index, data.iloc[:,1])]
        m = Pie()
        m.add('', plt_data)
        m.set_global_opts(legend_opts = opts.LegendOpts(is_show=False))
        m.set_series_opts(label_opts=opts.LabelOpts(is_show=False))
        return m
    
    def selected_static_describe(self, key):
        res = self.select[key].value_counts().to_frame()
        res['ratio'] = res/res.sum()
        return res
        
    def fast_group_counts(self, data, limits, key):
        if len(limits) < 1:
            res = data[key].value_counts()
            res = res.to_frame('num')
            res[key]=res.index
            res.index = range(res.shape[0])
            return res
        else:
            res = data.groupby(limits)[key].value_counts()    
        res = res.to_dict()

        k, v = list(res.keys()), list(res.values())

        k = [list(i) for i in k]

        res = []
        for i in zip(k,v):
            res.append(i[0]+[i[1]])
        if type(limits) != type([1,2,3]):
            res = pd.DataFrame(res, columns=[limits]+[key]+['num'])
        else:
            res = pd.DataFrame(res, columns=limits+[key]+['num'])
        return res
    
    def cal_exact_table(self, limit):
        if limit is None:
            limits = []
        else:
            limits = [limit]
        Kcol = self.fast_group_counts(self.data, limits+['yhaplogroup'], 'surname')

        Kcol.rename(columns = {'num':'K'} , inplace=True)

        Scol = self.fast_group_counts(self.data, limits, 'surname')
        Scol.rename(columns = {'num':'S'} , inplace=True)

        Ycol = self.fast_group_counts(self.data, limits, 'yhaplogroup')
        Ycol.rename(columns = {'num':'Y'} , inplace=True)

        result = pd.merge(Kcol,Ycol, on=limits+['yhaplogroup'], how='left')

        result = pd.merge(result, Scol, on=limits+['surname'], how='left')
        
        if limit is None:
            result['T'] = self.data.shape[0]
        else:
            Tcol = self.fast_group_counts(self.data, [], limit)
            Tcol.rename(columns = {'num':'T'}, inplace=True)
            result['T'] = pd.merge(result, Tcol, on=limit, how='left')['T']

        result['TT'] = result['K']
        result['TF'] = result['Y'] - result['K']
        result['FT'] = result['S'] - result['K']
        result['FF'] = (result['T']-result['Y'])- result['FT']
        result = result[limits+['yhaplogroup', 'surname', 'TT', 'TF', 'FT', 'FF']]
        valid = result[(result.TT != 0) &  (result.FT != 0) & (result.FF !=0)]  #(result.TF != 0) &
        tb = valid.iloc[:, 2:]
        f = lambda row: fisher_exact([[row['TT'], row['TF']], [row['FT'], row['FF']]], 'greater')[1]
        p = tb.apply(f, 1)
        valid['p_val'] = p
        return valid
    
    def add_prov_pop_density(self, tb, city_area_info):
        key = tb.keys().drop('province')[0]
        
        pop_info = city_area_info.groupby('province').sum()
        

        tb = pd.merge(tb, pop_info, on='province', how='left')

        tb['pop_density'] = tb['population']*10000*tb[key]/tb['area']

        tb.drop(['area', 'population'], axis=1, inplace=True)
        tb.rename(columns={'pop_density':'单倍群人口密度'}, inplace=True)
        return tb
    
    def add_city_pop_density(self, tmp, city_area_info):

        pop_info = city_area_info[['city', 'area']]

        tmp = pd.merge(tmp, pop_info, on='city', how='left')


        tmp['pop_density'] = tmp['hap_pop']*10000/tmp['area']

        tmp.drop(['area', 'population'], axis=1, inplace=True)
        tmp.rename(columns={'pop_density':'单倍群人口密度'}, inplace=True)
        return tmp

    
class TreeShape():
    def __init__(self, tree_data):
        self.tree = tree_data

    def parse_position(self, s):
        return re.findall('\d{4,}', s)
        
    def search_node_help(self, name, tree):
        if tree['name'] == name:
            return tree
        else:
            if tree['children'] is None:
                return None
            for i in tree['children']:
                target = self.search_node_help(name, i)
                if target is not None:
                    return target
                
    def search_node(self, hapname):
        return self.search_node_help(hapname, self.tree)
    
    def parse_pos(self, node):
        ids = node['rsids']
        pos = []
        for i in ids:
            if pd.isna(i):
                return []
            pos.append(re.findall('#(.*?):', i)[0])
        return pos
    
    def parse_pos_full(self, node):
        ids = node['rsids']
        pos = []
        for i in ids:
            if pd.isna(i):
                return []
            pos.append(re.findall('(.*?)#(.*?):(.*?)->(.*)', i)[0])
        pos = pd.DataFrame(pos, columns=['pos_name', 'pos', 'ref', 'alt'])
        pos['yhaplogroup'] = node['name']
        return pos
    
    def pull_pos_name_table(self, auto_save=True):
        res = []
        tree = self.tree
        def get_pos_name(tree):
            res.append(self.parse_pos_full(tree))
            for child in tree.get('children'):
                get_pos_name(child)
        get_pos_name(tree)
        res = pd.concat(res)
        if auto_save:
            res.to_csv('s3://genodig-prod/workspaces/users/ancestry/basic_data/pos_name_table.csv', index=False)
            print('pos_name_table.csv was auto saved at s3')
        return res
    
    def get_all_node_info(self, hapname):
        haps = self.get_branches(hapname)
        result = pd.DataFrame()
        for i in haps:
            result = result.append(self.parse_pos_full(self.search_node(i)))
        return result
    
    def get_pos(self, hapname):
        node = self.get_node(hapname)
        pos = self.parse_pos(node)
        return pos
    
    def search_pos_help(self, pos, tree):
        if pos in self.parse_pos(tree):
            return tree['name']
        else:
            if tree['children'] is None:
                return None
            for i in tree['children']:
                name = self.search_pos_help(pos, i)
                if name is not None:
                    return name
                
    def search_pos(self, pos):
        return self.search_pos_help(pos, self.tree)
                
    def search_barcode_help(self, barcode, tree):
        if 'barcode' in tree.keys():
            if tree['barcode'] == barcode:
                return tree
        else:
            for i in tree['children']:
                tmp = self.search_barcode_help(barcode, i)
                if tmp is not None:
                    return tmp
                
    def search_barcode(self, barcode):
        return self.search_barcode_help(barcode, self.tree)
    
    def get_parent_help(self, hapname, tree):
        if tree['children'] is None:
            return None
        for i in tree['children']:
            if i['name'] == hapname:
                return tree['name']
            else:
                tmp = self.get_parent_help(hapname, i)
                if tmp is not None:
                    return tmp
    
    def get_parent(self, hapname):
        return self.get_parent_help(hapname, self.tree)
    
    def get_path_help(self, hapname, tree):
        path = [hapname]
        parent = self.get_parent_help(hapname, tree)
        while parent is not None:
            path.append(parent)
            parent = self.get_parent_help(parent, tree)
        return path
    
    def get_path(self, hapname):
        return self.get_path_help(hapname, self.tree)
    
    def is_branch_of(self, son, parent):
        parent_node = self.search_node(parent)
        return self.search_node_help(son, parent_node) is not None
    
    def is_direct_branch_of(self, son, parent):
        dbranches = self.get_direct_branches(parent)
        return son in dbranches
    
    def get_branches(self, hapname):
        node = self.search_node(hapname)
        tmp = json.dumps(node)
        branches = re.findall('\"name\": \"(.*?)\"', tmp)
        return branches
    
    def get_direct_branches(self, hapname):
        node = self.search_node(hapname)
        dbranches = []
        for i in node['children']:
            dbranches.append(i['name'])
        return dbranches
    
    def get_time_help(self, hapname, tree):
        tmp = self.search_node_help(hapname, tree)
        parent = self.get_parent(hapname)
        if 'year' not in tmp.keys():
            return self.get_time_help(parent, tree)
        elif tmp['year'] == '':
            return self.get_time_help(parent, tree)
        else:
            return tmp['year']
    
    def get_time(self, hapname):
        return int(self.get_time_help(hapname, self.tree))
    
    def get_bc_time(self, hapname):
        t = int(self.get_time(hapname)) - 1950
        return t
    
    def cal_direct_branch_num(self, hapname):
        res = pd.DataFrame()
        branches = self.get_branches(hapname)
        res['yhaplogroup'] = branches

        time = []
        direct_branch_num = []
        for i in branches:
            time.append(self.get_time(i))
            direct_branch_num.append(len(self.get_direct_branches(i)))

        res['branch_num'] =direct_branch_num
        res['ybp'] = time
        res.sort_values(['ybp', 'branch_num'], ascending=False, inplace=True)
        return res
    
    def son_parent_pair_notification(self, df):
        tmp = df
        tmp.sort_values('ybp', inplace=True)
        son_parent_pair = []
        count = 0
        while count < (tmp.shape[0]-1):
            son = tmp.yhaplogroup.iloc[count]
            parent = tmp.yhaplogroup.iloc[count+1]
            if self.is_branch_of(son, parent):
                son_parent_pair.append('是 %s 的下游支系'%(parent))
            else:
                son_parent_pair.append('')
            count+=1
        son_parent_pair.append('')

        tmp['note'] = son_parent_pair

        tmp.sort_values('ybp', ascending=False, inplace=True)
        return tmp
    
    def son_parent_pair_notification2(self, df):
        tmp = df
        freq = tmp.ybp.value_counts()

        checks = freq[freq>1].keys()

        note = []
        for i in checks:
            tmp_df = tmp[tmp.ybp == i]
            found_list = []
            for hap1 in tmp_df.yhaplogroup.tolist():
                if hap1 in found_list: 
                    continue
                for hap2 in tmp_df.yhaplogroup.tolist():
                    if hap2 in found_list: continue
                    if self.is_direct_branch_of(hap1, hap2):
                        s = '是 %s 的下游支系'%(hap2)
                        found_list.append(hap1)
                        found_list.append(hap2)
                        note.append([hap1, s])
        note = pd.DataFrame(note,columns=['yhaplogroup', 'note'])

        new_tmp = fix_data(tmp, note, 'yhaplogroup', 'note')
        new_tmp.sort_values('ybp', ascending=False, inplace=True)
        return new_tmp
    
    def get_pos_name(self, pos):
        res = self.get_all_node_info(self.search_pos(pos))
        return res[res['pos'] == pos]
    
    def select_direct_branch_num_data(self,yhaplogroup, time_range, direct_branch_num_thresh, direct_branch_count_data):
        haps = self.get_branches(yhaplogroup)
        select = direct_branch_count_data[(direct_branch_count_data.yhaplogroup.isin(haps))
                                         &(direct_branch_count_data.ybp <= time_range[1])
                                         &(direct_branch_count_data.ybp >= time_range[0])
                                         &(direct_branch_count_data.branch_num >= direct_branch_num_thresh)]
        return select

    def visualize_direct_branch_num(self,hapname,time_range,thresh,branch_num_data=None):
        if branch_num_data is None:
            tmp = self.cal_direct_branch_num(hapname)
            tmp = self.select_direct_branch_num_data(hapname, time_range, thresh, branch_num_data)
        else:
            tmp = self.select_direct_branch_num_data(hapname, time_range, thresh, branch_num_data)
        plt.scatter(tmp.ybp, tmp.branch_num)
        top3 = tmp.sort_values('branch_num', ascending=False)
        top3 = top3.iloc[:3, :]
        used_branch_num = []
        for i in range(top3.shape[0]):
            text_info = top3.iloc[i, :]
            tmp_y = text_info['branch_num']
            while tmp_y in used_branch_num:
                tmp_y += 0.2
            used_branch_num.append(tmp_y)
            plt.text(text_info['ybp'], tmp_y, text_info['yhaplogroup'])
        plt.xlabel('ybp')
        plt.ylabel('direct_branch_num')
        return tmp
    
    def select_path(self, haplogroup, limit_ybp, branch_data):
        tmp = branch_data[branch_data.yhaplogroup.isin(self.get_path(haplogroup))]
        return tmp[tmp.ybp<=limit_ybp]
    
    def get_tmrc_node(self, haps):
        p0 = pd.DataFrame(self.get_path(haps[0]), columns=[haps[0]])
        order = p0.index.tolist()
        order.reverse()
        p0['order'] = order
        for i in haps[1:]:
            tmp = pd.DataFrame(self.get_path(i), columns=[i])
            order = tmp.index.tolist()
            order.reverse()
            tmp['order'] = order
            p0 = pd.merge(p0, tmp, on='order', how='outer')
        ks = p0.keys().drop('order')
        num = len(ks)
        name = p0[p0.apply(lambda x: x[ks].value_counts()[0]==num, axis=1)].sort_values('order', ascending=False)[haps[0]].iloc[0]
        n = self.search_node(name)
        return n
    
    def all_pos(self, pos, tree, container):
        if pos in self.parse_pos(tree):
            container.append(tree['name'])
        else:
            if tree['children'] is None:
                return None
            for i in tree['children']:
                self.all_pos(pos, i, container)

    def pos_freq(self, pos):
        res = []
        self.all_pos(pos, self.tree, res)
        return res
    
    def trash_pos_eval(self, position, node_name):
        res = []
        nodes = self.pos_freq(position)
        for i in nodes:
            if i != node_name:
                n = self.get_tmrc_node([node_name, i])
                res.append((i, n['year'], n['name']))
        res = pd.DataFrame(res, columns=['cmp_node', 'tmrca', 'ca_node'])
        return res
    
    def trash_node_eval(self,node_name):
        result = pd.DataFrame()
        for pos in self.parse_position(str(self.search_node(node_name)['rsids'])):
            tmp = self.trash_pos_eval(pos,node_name)
            tmp['pos'] = pos
            result = result.append(tmp)
        return result
    
    
    def get_end_nodes(self, node_name):
        def get_end_branches(node, container):
            if node['children'] == []:
                container.append(node['name'])
            else:
                for i in node['children']:
                    get_end_branches(i, container)
        container = []
        node = self.search_node(node_name)
        get_end_branches(node, container)
        return container

    def path_between(self, parent_node, son_node):
        tmp = self.get_path(son_node)
        if parent_node not in tmp:
            return None
        else:
            tmp = tmp[:tmp.index(parent_node)+1]
            tmp.reverse()
            return tmp
    
    def cmp_path(self, barcodes):

        result = []

        for i in barcodes:
            node = self.search_barcode(i)
            if node is None:
                print(i, 'not on ytree')
                continue
            name = node['name']
            path=self.get_path(name)
            path.reverse()
            num = len(path)
            result.append([name, path, num])

        result = pd.DataFrame(result, columns=['name', 'path', 'num'])

        valid_lenth = result.num.max()

        result = pd.DataFrame()
        for i in barcodes:
            node = self.search_barcode(i)
            if node is None:
                continue
            name = node['name']
            path=self.get_path(name)
            path.reverse()
            path = path + [None]*(valid_lenth-len(path))
            result[i]=path

        def tmr_node(row):
            keys = row.keys()
            tmp = row[0]
            for i in keys:
                if row[i] != row[0]:
                    return False
            return True

        result['is_common_node'] = result.apply(tmr_node, 1)
        display(result)
        print('最近公共节点:')
        print(result[result['is_common_node']].iloc[-1,0])
        return result
          
    def calculate_yfull_sample_result_table(self):
            sample_result_table = []
            def calculate_sample_result(tree):
                if 'barcode' in tree.keys():
                    barcode = tree['barcode']
                    name = tree['name']
                    sample_result_table.append([barcode, name])
                if 'children' in tree.keys():
                    for i in tree['children']:
                        calculate_sample_result(i)
            calculate_sample_result(self.tree)
            self.__sample_result = pd.DataFrame(sample_result_table, columns=['barcode', 'yhaplogroup'])
            
    @property
    def sample_result(self):
        try:
            return self.__sample_result
        except:
            self.calculate_yfull_sample_result_table()
            return self.__sample_result    
        
    def get_sub_samples(self, node_name):
        return self.sample_result[self.sample_result.yhaplogroup.isin(self.get_branches(node_name))]
    
    def search_pos_name(self, positions):
        pos_name = pd.DataFrame()
        for i in positions:
            pos_name = pos_name.append(self.get_pos_name(i))
        return pos_name
    
    def __get_chip_designed_haps(self):
        chip_designed_haps = dict()
        for i in self.all_data.version.unique():
            chip_designed_haps[i] = set(self.all_data[self.all_data.version == i].yhaplogroup.unique())
        self.__chip_designed_haps = chip_designed_haps
        
    @property
    def chip_designed_haps(self):
        try:
            return self.__chip_designed_haps
        except:
            self.__get_chip_designed_haps()
            return self.__chip_designed_haps
        
    def map_to_chip_version(self,hapname):    
        for version in self.chip_designed_haps.keys():
            def map_to_chip(haplogroup):
                if haplogroup is None:
                    return None
                if haplogroup in chip_designed_haps[version]:
                    return haplogroup
                else:
                    return map_to_chip(shape.get_parent(haplogroup))
            print(version,map_to_chip(hapname))
            
    def __cal_node_time(self):
        node_time = []     
        def get_node_time(tree, parent_time=1000000):
            name = tree.get('name')
            time = tree.get('year', parent_time)
            if time == '':
                time = parent_time
            time = int(time)
            node_time.append([name, time])
            for i in tree.get('children',[]):
                get_node_time(i, time)
        get_node_time(self.tree)
        node_time = pd.DataFrame(node_time, columns=['yhaplogroup', 'ybp'])
        self.__node_time = node_time
    
    @property
    def node_time(self):
        try:
            return self.__node_time
        except:
            self.__cal_node_time()
            return self.__node_time
        
    def generate_node_mapper(self, thresh):
        node_mapper = {}
        tree = self.tree
        def map_to_node_by_time(tree, thresh=2300, map_node='ROOT', parent_time=10000):
            name, time = tree.get('name'), tree.get('year', parent_time)
            if time=='':
                time=parent_time
            time = int(time)
            if (parent_time>thresh) & (time<=thresh):
                for i in tree.get('children'):
                    map_to_node_by_time(i, thresh, name, time)
            elif parent_time>thresh:
                for i in tree.get('children'):
                    map_to_node_by_time(i, thresh, name, time)
            else:
                node_mapper[name] = map_node
                for i in tree.get('children'):
                    map_to_node_by_time(i, thresh, map_node, time)
        map_to_node_by_time(tree, thresh=thresh)
        return node_mapper
    
    def get_full_path(self, yhaplogroup):
        return self.path_between('ROOT', yhaplogroup)[:-1]+self.get_branches(yhaplogroup)
    
    
# csv数据库

class CsvDataBase():
    def create_db(self, cols, dtypes):
        self.db = pd.DataFrame(columns=cols)
        self.dtypes = dtypes

    def save(self, filename):
        for i in range(len(self.dtypes)):
            if self.dtypes[i] == 'l':
                key = self.db.keys()[i]
                self.db[key] = self.db[key].map(lambda x: self.list_to_str(x))
        self.db.loc['dtypes'] = self.dtypes
        self.db.to_csv(filename + '.csv')
        print('saving complete.')

    def load(self, filename):
        self.db = pd.read_csv(filename + '.csv', index_col=0)
        self.dtypes = self.db.loc['dtypes'].tolist()
        self.delete_by_loc('dtypes')
        for i in range(len(self.dtypes)):
            if self.dtypes[i] == 'l':
                key = self.db.keys()[i]
                self.db[key] = self.db[key].map(lambda x: self.str_to_list(x))
            else:
                try:
                    key = self.db.keys()[i]
                    self.db[key] = self.db[key].map(lambda x: eval(x))
                except:
                    pass
        print('loading complete.')

    def list_to_str(self, L):
        if L is None:
            return 'None'
        s = ''
        for i in L:
            s += str(i) + ','
        return s[:-1]

    def str_to_list(self, S):
        if S == 'None':
            return None
        tmp = S.split(',')
        try:
            tmp = [eval[i] for i in tmp]
        except:
            pass
        return tmp

    def add_item(self, vals):
        self.db.loc[self.db.shape[0] + 1] = vals

    def delete_by_loc(self, index):
        self.db.drop(index, inplace=True)

    def delete_by_condition(self, condition):
        indexs = self.db[condition].index.tolist()
        for i in indexs:
            self.delete_by_loc(i)

    def update(self, index, vals):
        self.delete_by_loc(index)
        self.db.loc[index] = vals
        
        
class SimpleData():
    def __init__(self, data_frame, tree_shape):
        '''
        data_frame: pandas.DataFrame
        tree_shape: Zero defined TreeShape object
        '''
        self.all_data = data_frame #all_data
        self.shape = tree_shape

    
    def initialize_data_set(self):
        self.data = self.all_data.copy(deep=True) # data set used to calculate
        print('data set size: %s'%(self.data.shape[0]))
    
    def initialize_selected_data(self):
        # seleceted data for calculation
        self.selected_data = self.all_data.copy(deep=True)
        print('selected data size: %s'%(self.selected_data.shape[0]))
        
    def data_set_selector(self, **kwargs):
        for key in kwargs.keys():
            val = kwargs[key]
            if type(val) == type([1]):
                self.data = self.data[self.data[key].isin(val)]
            else:
                self.data = self.data[self.data[key] == val]
        print('data set size: %s'%(self.data.shape[0]))
        
        
    def change_hap_name_in_chip_version(self, chip_version, chip_hap_name, to_new_hap_name):
        self.data['yhaplogroup'][(self.data.version==chip_version) & (self.data.yhaplogroup == chip_hap_name)] = to_new_hap_name
        print('change %s to %s in version %s'%(chip_hap_name, to_new_hap_name, chip_version))
    
    
    def data_selector(self, **kwargs):
        for key in kwargs.keys():
            val = kwargs[key]
            if type(val) == type([1]):
                self.selected_data = self.selected_data[self.selected_data[key].isin(val)]
            else:
                self.selected_data = self.selected_data[self.selected_data[key] == val]
#         self.data_set_selector(version=self.selected_data.version.unique().tolist())
#         print('version: ', self.selected_data.version.unique().tolist(), 'are used')
        print('selected data size: %s'%(self.selected_data.shape[0]))
    

    def get_data(self, **kwargs):
        tmp = self.data
        for key in kwargs.keys():
            val = kwargs[key]
            if type(val) == type([1]):
                tmp = tmp[tmp[key].isin(val)]
            else:
                tmp = tmp[tmp[key] == val]
        return tmp
    
    
    def get_data_from_all_data(self, **kwargs):
        tmp = self.all_data
        for key in kwargs.keys():
            val = kwargs[key]
            if type(val) == type([1]):
                tmp = tmp[tmp[key].isin(val)]
            else:
                tmp = tmp[tmp[key] == val]
        return tmp

    def calculate_yhaplogroup_surname_exact_test_table(self, yhaplogroup, surname):
        yhaplogroup_surname = self.data[(self.data.yhaplogroup == yhaplogroup) & (self.data.surname == surname)].shape[0]
        yhaplogroup_Nsurname = self.data[(self.data.yhaplogroup == yhaplogroup) & (self.data.surname != surname)].shape[0]
        Nyhaplogroup_surname = self.data[(self.data.yhaplogroup != yhaplogroup) & (self.data.surname == surname)].shape[0]
        Nyhaplogroup_Nsurname = self.data[(self.data.yhaplogroup != yhaplogroup) & (self.data.surname != surname)].shape[0] 
        return [[yhaplogroup_surname, yhaplogroup_Nsurname], [Nyhaplogroup_surname, Nyhaplogroup_Nsurname]]
    
    def calculate_yhaplogroup_surname_exact_test_table_include_sub_branches(self, yhaplogroup, surname):
        yhaplogroup_surname = self.data[(self.data.yhaplogroup.isin(self.shape.get_branches(yhaplogroup))) & (self.data.surname == surname)].shape[0]
        yhaplogroup_Nsurname = self.data[(self.data.yhaplogroup.isin(self.shape.get_branches(yhaplogroup))) & (self.data.surname != surname)].shape[0]
        Nyhaplogroup_surname = self.data[(self.data.yhaplogroup.isin(self.shape.get_branches(yhaplogroup))==False) & (self.data.surname == surname)].shape[0]
        Nyhaplogroup_Nsurname = self.data[(self.data.yhaplogroup.isin(self.shape.get_branches(yhaplogroup))==False) & (self.data.surname != surname)].shape[0] 
        return [[yhaplogroup_surname, yhaplogroup_Nsurname], [Nyhaplogroup_surname, Nyhaplogroup_Nsurname]]
    
    
    
    def calculate_fisher_exact_p_value(self, exact_table):
        '''
        exact_table:
            [
                [yhaplogroup_surname, yhaplogroup_Nsurname], 
                [Nyhaplogroup_surname, Nyhaplogroup_Nsurname]
            ]
        like table
        '''
        return fisher_exact(exact_table, 'greater')[1]
    
    def selected_data_ratio_over_data_set(self, key):
        '''
            key: province, city, county
            calculate selected_data/data_set on key
        '''
        use_version = self.selected_data.version.unique().tolist()
#         self.data_set_selector(version=use_version)
        sample_freq = self.selected_data[key].value_counts()
        sample_df = pd.DataFrame()
        sample_df[key] = sample_freq.keys().tolist()
        sample_df['sample_num'] = sample_freq.values
#         display(sample_df)
        total_freq = self.data[self.data.isin(sample_freq.keys().tolist())][key].value_counts()
        total_df = pd.DataFrame()
        total_df[key] = total_freq.keys().tolist()
        total_df['total_num'] = total_freq.values
#         display(total_df)
        ratio_table = pd.merge(sample_df, total_df, on=key, how='left')
        ratio_table['ratio'] = ratio_table['sample_num']*1.0/ratio_table['total_num']
        ratio_table.sort_values('ratio', ascending=False, inplace=True)
        return ratio_table

    def population_ratio_in_china(self, key_population_weight_info, weight_by='city', return_full_table=False):
        '''
            加权计算选中数据在中国人口占比，默认按city加权
        '''
        ratio_table = self.selected_data_ratio_over_data_set(weight_by)
        if weight_by == 'city':
#             weight_info = key_population_weight_info[key_population_weight_info[weight_by].isin(ratio_table[weight_by])]#[[weight_by, 'weight']]
            ratio_table = pd.merge(ratio_table, key_population_weight_info, on=weight_by, how='outer')
            ratio_table['ratio'] = ratio_table.sample_num/ratio_table.total_num
            ratio_table['hap_pop'] = ratio_table.ratio*ratio_table.population
#             ratio_table['weighted_ratio'] = ratio_table.hap_pop/ratio_table.population
            if return_full_table:
                return ratio_table
            return ratio_table.hap_pop.sum()/ratio_table.population.sum()
        else:
            total_pop = key_population_weight_info.population.sum()
            weight_info = key_population_weight_info.groupby(weight_by)['population'].sum()
            weight_df = pd.DataFrame()
            weight_df['province'] = weight_info.index.tolist()
            weight_df['population'] = weight_info.values
            weight_df['weight'] = weight_df['population']/total_pop
            ratio_table = pd.merge(ratio_table, weight_df, on=weight_by, how='outer')
            ratio_table['ratio'] = ratio_table.sample_num/ratio_table.total_num
            ratio_table['hap_pop'] = ratio_table.ratio*ratio_table.population
#             ratio_table['weighted_ratio'] = ratio_table.hap_pop/ratio_table.population
            if return_full_table:
                return ratio_table
            return ratio_table.hap_pop.sum()/ratio_table.population.sum()


    def draw_province_distribution_map(self, weight, weight_by='province', weighted_ratio=True):
        tb = self.population_ratio_in_china(weight, weight_by, return_full_table=True)
        tb = tb.groupby('province')[['sample_num', 'total_num', 'hap_pop', 'population']].sum()
        tb['ratio'] = tb.sample_num/tb.total_num
        tb['weighted_ratio'] = tb.hap_pop/tb.population
        tb.sort_values('weighted_ratio', ascending=False, inplace=True)
        tb.fillna(0,inplace=True)
        tb['province'] = tb.index
        tb.index = range(tb.shape[0])
        if weighted_ratio:
            tb = tb[['province', 'weighted_ratio']]
            tb.rename(columns={'weighted_ratio':'ratio'},inplace=True)
        tb = tb[['province','ratio']]
        max_val, min_val = tb.iloc[:, 1].max(), tb.iloc[:, 1].min()
        map_data = []
        for i in range(tb.shape[0]):
            tmp = tb.iloc[i, :]
            map_data.append([tmp.iloc[0], tmp.iloc[1]])
        DIC = {True:'（含下游）', False:'（不含下游）'}
        m = Map()
        m.add('', map_data, "china", is_map_symbol_show=False, is_selected=False)
        m.set_global_opts(
            toolbox_opts=opts.ToolBoxFeatureOpts(save_as_image='png'),
            title_opts=opts.TitleOpts(title='分布频率'),
            visualmap_opts=opts.VisualMapOpts(
                is_show=False,
                orient='horizontal',\
                max_=max_val, min_=min_val, \
                range_color=['#FFFFFF','#C4DCF9','#99CAFF','#66B0FF','#3395FF'],\
                range_text = ['占比较大', '占比较小'],
#                 split_number=5,\
#                 is_piecewise=True
            ),
            legend_opts=opts.LegendOpts(is_show=False)
        )
        m.set_series_opts(
            label_opts=opts.LabelOpts(is_show=False)
        )
        tb.sort_values('ratio', ascending=False, inplace=True)
        return m, tb
    

    def add_prov_pop_density(self, tb, city_area_info):
        key = tb.keys().drop('province')[0]


        pop_info = city_area_info.groupby('province').sum()


        tb = pd.merge(tb, pop_info, on='province', how='left')

        tb['pop_density'] = tb['population']*10000*tb[key]/tb['area']

        tb.drop(['area', 'population'], axis=1, inplace=True)
        tb.rename(columns={'pop_density':'单倍群人口密度'}, inplace=True)
        return tb
    
    
    def add_city_pop_density(self, tmp, city_area_info):

        pop_info = city_area_info[['city', 'area']]

        tmp = pd.merge(tmp, pop_info, on='city', how='left')


        tmp['pop_density'] = tmp['hap_pop']*10000/tmp['area']

        tmp.drop(['area', 'population'], axis=1, inplace=True)
        tmp.rename(columns={'pop_density':'单倍群人口密度'}, inplace=True)
        return tmp

    
    def draw_map_from_input_file(self, tb, weight):
        tb = pd.merge(tb, weight['province'].drop_duplicates(), on='province', how='outer')
        tb.fillna(0,inplace=True)
        tb = tb[['province', 'ratio']]
        max_val, min_val = tb.iloc[:, 1].max(), tb.iloc[:, 1].min()
        map_data = []
        for i in range(tb.shape[0]):
            tmp = tb.iloc[i, :]
            map_data.append([tmp.iloc[0], tmp.iloc[1]])
        DIC = {True:'（含下游）', False:'（不含下游）'}
        m = Map()
        m.add('', map_data, "china", is_map_symbol_show=False, is_selected=False)
        m.set_global_opts(
            toolbox_opts=opts.ToolBoxFeatureOpts(save_as_image='png'),
            title_opts=opts.TitleOpts(title='分布频率'),
            visualmap_opts=opts.VisualMapOpts(
                is_show=False,
                orient='horizontal',\
                max_=max_val, min_=min_val, \
                range_color=['#FFFFFF','#C4DCF9','#99CAFF','#66B0FF','#3395FF'],\
                range_text = ['占比较大', '占比较小'],
#                 split_number=5,\
#                 is_piecewise=True
            ),
            legend_opts=opts.LegendOpts(is_show=False)
        )
        m.set_series_opts(
            label_opts=opts.LabelOpts(is_show=False)
        )
        return m
    
    
    def draw_city_distribution_map(self, weight, province, weighted_ratio=True):
        numerator = self.selected_data[self.selected_data.province==province]
        denominator = self.data[self.data.province==province]
        numerator = numerator.city.value_counts().to_frame('num')
        denominator = denominator.city.value_counts().to_frame('total_num')
        tb = pd.concat((numerator, denominator), axis=1)
        tb['ratio'] = tb.num/tb.total_num
        tb['city'] = tb.index
        tb.index=range(tb.shape[0])
        tb.fillna(0, inplace=True)
        tb = tb[['city', 'num', 'total_num', 'ratio']]
        weight = weight[weight.province==province]
        tb = pd.merge(tb, weight, on='city', how='left')
        tb['weighted_ratio'] = tb.ratio*tb.weight
        if weighted_ratio:
            tb = tb[['city', 'weighted_ratio']]
            tb.rename(columns={'weighted_ratio':'ratio'},inplace=True)
        tb = tb[['city','ratio']]
        max_val, min_val = tb.iloc[:, 1].max(), tb.iloc[:, 1].min()
        map_data = []
        for i in range(tb.shape[0]):
            tmp = tb.iloc[i, :]
            map_data.append([tmp.iloc[0], tmp.iloc[1]])
        DIC = {True:'（含下游）', False:'（不含下游）'}
        m = Map()
        m.add('', map_data, province, is_map_symbol_show=False, is_selected=False)
        m.set_global_opts(
            title_opts=opts.TitleOpts(title='分布频率'),
            visualmap_opts=opts.VisualMapOpts(
                is_show=False,
                orient='horizontal',\
                max_=max_val, min_=min_val, \
                range_color=['#FFFFFF','#C4DCF9','#99CAFF','#66B0FF','#3395FF'],\
                range_text = ['占比较大', '占比较小'],
#                 split_number=5,\
#                 is_piecewise=True
            ),
            legend_opts=opts.LegendOpts(is_show=False)
        )
        m.set_series_opts(
            label_opts=opts.LabelOpts(is_show=False)
        )
        return m, tb
    
    def draw_composition_pie_chart(self):
        res = self.selected_data[composition_key].sum()
        res = res/self.selected_data.shape[0]
        res = res[res>0]
        plt_data = [list(i) for i in zip(res.keys(), res.values)]
        m = Pie()
        m.add('', plt_data)
        m.set_global_opts(title_opts=opts.TitleOpts(title="民族成分"), legend_opts = opts.LegendOpts(is_show=False))
        return (m,res)
        
    def __get_yhaplogroup_info(self, barcodes):
        self.shape.calculate_yfull_sample_result_table()
        ysample = self.shape.sample_result[self.shape.sample_result.barcode.isin(barcodes)]
        
        print('yfull result:')
        print(list(ysample.barcode))
        chip_sample = self.get_data_from_all_data(barcode=barcodes)[['barcode', 'yhaplogroup']]
        chip_sample = chip_sample[(chip_sample.yhaplogroup.notna()) & (chip_sample.yhaplogroup!='') & (chip_sample.yhaplogroup!='female')]
        
        print('chip result:')
        print(chip_sample.barcode.tolist())
        valid_hap = ysample.append(chip_sample)
        
        print('not found:')
        not_found_bars = set(barcodes).difference(set(valid_hap.barcode))
        print(not_found_bars)
        valid_hap.drop_duplicates('barcode',keep='first', inplace=True)
        valid_hap.sort_values('yhaplogroup', inplace=True)
        return valid_hap
    
    def __generate_cmp_path_table(self, hap_list):
        hap_list = list(set(hap_list))
        table = pd.DataFrame()
        table['name'] = hap_list
        table['yhaplogroup'] = hap_list
        return table
    
    def cmp_path(self, barcodes=None, hap_list=None, mix_query=False):
        if mix_query:
            if barcodes is not None:
                valid_hap = self.__get_yhaplogroup_info(barcodes)
            if hap_list is not None:
                valid_hap = self.__generate_cmp_path_table(hap_list)
        else:
            valid_hap = self.shape.sample_result[self.shape.sample_result.barcode.isin(barcodes)]
        result = []

        for name in valid_hap.yhaplogroup.tolist():
            path=self.shape.get_path(name)
            path.reverse()
            num = len(path)
            result.append([name, path, num])

        result = pd.DataFrame(result, columns=['name', 'path', 'num'])


        valid_lenth = result.num.max()
        result1=pd.DataFrame()
        for i in range(result.shape[0]):
            barcode = valid_hap.iloc[i,0]
            name = result.iloc[i,0]
            path = result.iloc[i,1]
            path = path + [None]*(valid_lenth-len(path))
            result1[barcode]=path

        result1['common_num'] = result1.apply(lambda x: x.value_counts()[0], 1)
        display(result1)
        print('最近公共节点:')
        tmrcn = result1[result1['common_num'] == result1['common_num'].max()].iloc[-1, 0]
        print(tmrcn)
        print('最近共祖时间: %s ybp'%(self.shape.search_node(tmrcn).get('year', 0)))
#         return result1
        return tmrcn
        
    def exclude_values(self, key, values):
        tmp=self.all_data[key]
        if type(values) == type(list()):
            return list(set(tmp[tmp.isin(values)==False]))
        else:
            return list(set(tmp[tmp!=values]))
        
    def include_sub_haps(self, yhaplogroup):
        return self.shape.get_branches(yhaplogroup)
    
    def get_tmrc_non_empty_node(self, yhaplogroup, province, surname):
        iter_yhap = yhaplogroup
        bars = self.shape.get_sub_samples(iter_yhap).barcode.tolist()
        tmp = self.get_data_from_all_data(barcode=bars, surname=surname, province=province)
        year = self.shape.search_node(iter_yhap).get('year', 0)
        while tmp.shape[0] < 1:
            iter_yhap = self.shape.get_parent(iter_yhap)
            bars = self.shape.get_sub_samples(iter_yhap).barcode.tolist()

            tmp = self.get_data_from_all_data(barcode=bars, surname=surname, province=province)
            year = self.shape.search_node(iter_yhap).get('year', 0)
        #     print(iter_yhap, year, tmp.shape[0])
            if iter_yhap == 'ROOT':
                break

        if iter_yhap != 'ROOT':
            print(iter_yhap, year)
            return tmp

    def generate_exact_matrix(self, key1, key2):
        '''
        对任意两个特征进行显著性检验
        '''
        try:
            self.selected_data
        except:
            print('未选择数据')
            return None

        key1_key2_gp = self.selected_data.groupby([key1, key2]).barcode.count()

        result = pd.DataFrame(key1_key2_gp.index.tolist(), columns=[key1, key2])
        

        result['TT'] = key1_key2_gp.values.tolist()

        key1_df = self.selected_data[key1].value_counts().to_frame('key1_count')
        key1_df.index.name = key1

        key2_df = self.selected_data[key2].value_counts().to_frame('key2_count')
        key2_df.index.name = key2

        result = pd.merge(result, key1_df, on=key1, how='left')
        result = pd.merge(result, key2_df, on=key2, how='left')

        result['total'] = self.selected_data.shape[0]

        result['FT'] = result.key1_count - result.TT
        result['TF'] = result.key2_count - result.TT
        result['FF'] = result.total - result.key1_count - result.key2_count + result.TT

        result[[key1, key2, 'TT', 'FT', 'TF', 'FF']]

        f = lambda row: fisher_exact([[row['TT'], row['FT']], [row['TF'], row['FF']]], 'greater')[1]
        result['p_val'] =  result.apply(f, 1)
        # 筛选出p<0.05的数据，认为这些姓氏与单倍群具有显著性关联。阈值可调到0.01结果更可靠但覆盖用户会减少
        exact_result = result[(result['p_val'] < 0.05) & (result.TT>=3)].sort_values('p_val', ascending=True)

        return exact_result[[key1, key2, 'TT', 'FT', 'TF', 'FF', 'p_val']]
    

class ChineseNumber():    
    def __init__(self, number):
        self.num = number
        self.num_to_chinese = {
            0:'零',
            1:'一',
            2:'两',
            3:'三',
            4:'四',
            5:'五',
            6:'六',
            7:'七',
            8:'八',
            9:'九',
        }

        self.digits_to_chinese = {
                10:'十',
                100:'百',
                1000:'千',
                10000:'万',
                100000:'十万',
                1000000:'百万',
                10000000:'千万',
        }
    
    def parse_number(self):
        if self.num < 1:
            for i in self.digits_to_chinese:
                if (self.num * i <= 10) & (self.num * i > 1):
                    self.real_part = round(self.num * i,1)
                    self.digits = i
                    self.positive=False
                    break
        else:
            for i in self.digits_to_chinese:
                if (self.num / i < 10) & (self.num / i >= 1):
                    self.real_part = round(self.num / i, 1)
                    self.digits = i
                    self.positive=True
                    break

    def describe_number(self):
        self.parse_number()
        if self.positive:
            _1st = str(self.real_part)[0]
            _2nd = str(self.real_part)[-1]
            s = self.num_to_chinese[int(_1st)] + self.digits_to_chinese[self.digits] + self.num_to_chinese[int(_2nd)] + self.digits_to_chinese[self.digits/10]
        else:
            denominator = int(round(self.digits/self.real_part, 0))
            _1st = str(denominator)[0]
            _2nd = str(denominator)[1]
            s = self.num_to_chinese[int(_1st)] + self.digits_to_chinese[self.digits/10] + self.num_to_chinese[int(_2nd)] + self.digits_to_chinese[self.digits/100] + '分之一  ' + str(denominator)
        return s
    
    def __repr__(self):
        s = str(self.num)+'\t'
        s += self.describe_number()
        return s

class PosQuery():
    def __init__(self, positions=[]):
        self.positions = positions
        self.queried_pos = set()
        self.query_raw_result = pd.DataFrame(columns=['s', 'locus.contig', 'locus.position', 'alleles', 'genotype', 'uid',
       'version'])
    
    def add_position(self, new_positions):
        '''添加位点到查询列表'''
        self.positions += new_positions
        
    def threading_pos_query(self):
        '''开新线程查询位点突变情况'''
        t = threading.Thread(target=self.pos_query, name='pos_query')
        t.start()
        if not t.isAlive:
            t.join()
            print(self.queried_pos, 'done')
    
    def pos_query(self):
        pos = list(set(self.positions).difference(self.queried_pos))
        self.queried_pos = self.queried_pos.union(set(self.positions).difference(self.queried_pos))
        if len(pos) < 1:
            return None
        print('query position:\n', pos)
        pos = ['Y:'+str(i) for i in pos]
        mt = (
            geno('affy')
            .locus(*pos)
            .filter(
                fields.gt.computed_gender == "male",
                fields.gt.version.inside('affy-v2.2', 'affy-v2.1', 'affy-2', 'affy-1'), 
                fields.gt.GT==hl.call(1,1)
            )
            .data())
        ht_c = mt.cols().flatten().select("s", "uid", "version").key_by('s')
        ht_e = mt.select_entries("genotype").key_cols_by().entries().select("s", "genotype").key_by('s')
        ht = ht_e.join(ht_c, how='left')
        df = ht.to_pandas()
        df = df[df.uid != '-1']
        df = df[df.genotype != './.']

        mt = (
            geno('ilmn')
            .locus(*pos)
            .filter(
                fields.ilmn_gt.gender == "M",
    #             fields.gt.version.inside('affy-v2.2', 'affy-v2.1', 'affy-2', 'affy-1'), 
                fields.gt.GT==hl.call(1,1)
            )
            .data())
        ht_c = mt.cols().flatten().select("s", "uid", "version").key_by('s')
        ht_e = mt.select_entries("genotype").key_cols_by().entries().select("s", "genotype").key_by('s')
        ht = ht_e.join(ht_c, how='left')
        df1 = ht.to_pandas()
        df1 = df1[df1.uid != '-1']
        df1 = df1[df1.genotype != './.']

        df = pd.concat([df, df1])
        self.query_raw_result = pd.concat([self.query_raw_result, df])
    
    @property
    def samples_count(self):
        '''统计每个位点在每版芯片的突变样本数量情况'''
        self.formalize_query_result()
        def uid_count(s):
            if s is None:
                return 0
            else:
                return len(s)
        tmp = self.query
        tmp[tmp.keys()[1:]] = tmp[tmp.keys()[1:]].applymap(uid_count)
        return tmp
    
    @property
    def query_result(self):
        '''统计每个位点在每版芯片有哪些位点突变'''
        self.formalize_query_result()
        return self.query
    
    def threading_is_all_query_done(self):
        if 'finish check' in str(threading.enumerate()):
            return None
        '''检查是否所有查询已经完成'''
        def f():
            s = time.time()
            while not self.is_all_query_done():
                e = time.time()
                if e-s > 3000:
                    print('timeout')
                    break
            print('ALL QUERY DONE')
        t = threading.Thread(target=f, name='finish check')
        t.start()
        
    def is_all_query_done(self):
        return 'pos_query' not in str(threading.enumerate())
    
    def formalize_query_result(self):
        res = self.query_raw_result
        query = []
        for p in res['locus.position'].unique().tolist():
            tmp_df = res[res['locus.position'] == int(p)]
            tmp_uid = tmp_df.uid.tolist()
            v1 = tmp_df[tmp_df['version'] == 'affy-1']
            v1 = v1.uid.tolist()
            if len(v1) < 1:
                v1 = None
            v2 = tmp_df[tmp_df['version'] == 'affy-2']
            v2 = v2.uid.tolist()
            if len(v2) < 1:
                v2 = None
            v2_1 = tmp_df[tmp_df['version'] == 'affy-v2.1']
            v2_1 = v2_1.uid.tolist()
            if len(v2_1) < 1:
                v2_1 = None
            v2_2 = tmp_df[tmp_df['version'] == 'affy-v2.2']
            v2_2 = v2_2.uid.tolist()
            if len(v2_2) < 1:
                v2_2 = None

            ilmn1_0 = tmp_df[tmp_df['version'] == 'ilmn-v1.0']
            ilmn1_0 = ilmn1_0.uid.tolist()
            if len(ilmn1_0) < 1:
                ilmn1_0 = None

            ilmn1_1 = tmp_df[tmp_df['version'] == 'ilmn-v1.1']
            ilmn1_1 = ilmn1_1.uid.tolist()
            if len(ilmn1_1) < 1:
                ilmn1_1 = None
            query.append([p, v1, v2, v2_1, v2_2, ilmn1_0, ilmn1_1])
        query = pd.DataFrame(query, columns = ['pos', 'affy-1', 'affy-2', 'affy-v2.1', 'affy-v2.2', 'ilmn-v1.0', 'ilmn-v1.1'])
        query['pos'] = query['pos'].astype(str)
        self.query = query
        
    def get_mut_samples(self, pos, version):
        pos = int(pos)
        version_map = {
            'v1' : 'affy-1',
            'v2' : 'affy-2', 
            'v2.1' : 'affy-v2.1',
            'v2.2' : 'affy-v2.2',
            'i1' : 'ilmn-v1.0',
            'i1.1' : 'ilmn-v1.1'
        }
        version = [version_map.get(i) for i in version]
        uids = set(self.query_raw_result[(self.query_raw_result['locus.position']==pos) & (self.query_raw_result.version.isin(version))].uid)
        barcode = [uid.uid_int_to_barcode(i) for i in uids]
        return barcode
    
    def cmp_pos(self, pos1, pos2, version):
        s1 = set(self.get_mut_samples(pos1, version))
        s2 = set(self.get_mut_samples(pos2, version))
        inner = len(s1.intersection(s2))
        if inner < 1:
            print('intersection ratio: %s' % (0))
        else:
            print('intersection ratio: %s' % (inner/len(s1)))
            
    def add_pos_name_info(self, pos_name_table):
        res = pd.merge(self.samples_count, pos_name_table, on='pos', how='left')
        return res
            
            
class DataLoader:
    def __init__(self):
        pass
    
    def load_yfull_tree(self):
        '''
        加载全序树形
        :return: tree, handler, haploTree
        '''
        from genodig.analysis.haploGroup import HaplogroupTree
        from genodig.analysis.haplo_edit import MergeTree
        from genodig.analysis.haplo_chip import MappingChip
        from genodig.analysis.halplo_years import CalcYears
        from genodig.analysis.haplo_single import SingleSample_Ontree
        from genodig.common.timeutil import time_now_str
        # 加载模块，初始化数据
        handler = SingleSample_Ontree()

        haploTree = HaplogroupTree(tree_object=handler.full_tree)

        # 加载单倍群树数据
        tree = haploTree.tree[0]
        return tree, handler, haploTree
    
    def load_pop_weight(self):
        '''
        加载省份人口权重数据
        '''
        return pd.read_csv('s3://genodig-prod/workspaces/users/ancestry/basic_data/city_population.csv')
    
    def load_male_data(self):
        '''加载清洗后男性用户数据'''
        return pd.read_csv('s3://genodig-prod/workspaces/users/ancestry/basic_data/male_data.csv')
    
    def load_all_data(self):
        '''加载全部用户数据'''
        return pd.read_csv('s3://genodig-prod/workspaces/users/ancestry/basic_data/all_data.csv')
    
    def load_pos_name_table(self):
        '''加载位点名称以及所在全序树节点名称'''
        return pd.read_csv('s3://genodig-prod/workspaces/users/ancestry/basic_data/pos_name_table.csv')
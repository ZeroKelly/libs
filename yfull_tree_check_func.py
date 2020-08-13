import re
import json
import pandas as pd
from genodig.common import uid
from genodig.common import blackhole
from genodig.common.fs import latest_path

import matplotlib.pyplot as plt

bar_pattern = re.compile('\d{3}-\d{4}-\d{4}')
qc_pattern = re.compile('\((.*?)[\|,](.*?)\)')


def hap_certainty(barcodes, yloci):
    loci = [ "	Y	"+ylocus for ylocus in yloci]
    uids = [ uid.barcode_to_uid(barcode) for barcode in barcodes]
    mapping = {}
    for i in uids:
        mapping[i] = uid.uid_int_to_barcode(i)
    gt_df = blackhole.get_table_data("yfullGt", uids, loci)
    gt_df.rename(columns = mapping, inplace=True)
    qt_df = blackhole.get_table_data("yfullQt", uids, loci)
    qt_df.rename(columns = mapping, inplace=True)
    mix_df = blackhole.yfull_raw_quality(uids, yloci)
    mix_df.rename(columns = mapping, inplace=True)
    res = mix_df
    missing=set(barcodes).difference(set(res.keys()[2:])) 
    for i in missing:
        res[i] = 'None'
    res = res.reindex(yloci)
    return res[barcodes]

def parse_barcode(s):
    return re.findall('\d{3}-\d{4}-\d{4}', s)

def parse_yhaplogroup(s):
    return re.findall('([A-Z]+?\-[A-Z]+\d+)', s)

def parse_surname(s):
    return re.findall('(.)氏', s)

def parse_position(s):
    return re.findall('\d{4,}', s)

def parse_rsids(node):
    if node is None:
        return None
    res = []
    for i in node['rsids']:
        pos, ref, alt = re.findall('#(.*?):(.*?)->(.*)', i)[0]
        res.append((pos, ref, alt))
    return res

def search_name(name, tree):
    if tree['name'] == name:
        return tree
    else:
        if tree['children'] is None:
            return None
        for i in tree['children']:
            target = search_name(name, i)
            if target is not None:
                return target
            
def search_barcode(barcode, tree):
    if 'barcode' in tree.keys():
        if tree['barcode'] == barcode:
            return tree
    else:
        if tree['children'] is None:
            return None
        for i in tree['children']:
            target = search_barcode(barcode, i)
            if target is not None:
                return target
            
def search_node(name, tree):
    if len(parse_barcode(name)) > 0:
        return search_barcode(name, tree)
    else:
        return search_name(name, tree)

def pos_comparison(pos_list1, pos_list2):
    s1 = set(pos_list1)
    s2 = set(pos_list2)
    if s1 != s2:
        redundent = s1.difference(s2)
        missing = s2.difference(s1)
#         print('自动树缺失位点:\n')
#         for i in redundent:
#             print(i)
#         print('\n自动树多出位点:\n')
#         for i in missing:
#             print(i)
        return (list(redundent), list(missing))
    return None,None

def generate_report(pos1, pos2):
    redundent, missing = pos_comparison(pos1, pos2)
    r = pd.DataFrame(redundent, columns=['position'])
    r['自动树状态'] = '缺少'
    m = pd.DataFrame(missing, columns=['position'])
    m['自动树状态'] = '多余'
    r = r.append(m)
    tmp = r.position.tolist()
    pos_on_tree = []
    pos_on_tree1 = []
    for i in tmp:
        pos_on_tree.append(search_pos(i,tree))
        pos_on_tree1.append(search_pos(i, manual_tree))
    r['在自动树上节点位置'] = pos_on_tree
    r['在人工树上节点位置'] = pos_on_tree1
    return r

def parse_pos(node):
    ids = node['rsids']
    pos = []
    for i in ids:
        if pd.isna(i):
            return []
        pos.append(re.findall('#(.*?):', i)[0])
    return pos

def search_pos(pos, tree):
    if pos in parse_pos(tree):
        return tree['name']
    else:
        if tree['children'] is None:
            return None
        for i in tree['children']:
            name = search_pos(pos, i)
            if name is not None:
                return name


def flatten_tree(container, tree):
    if tree['children'] is None:
        return None
    else:
        container.append(tree)
        for i in tree['children']:
            tmp = flatten_tree(container, i)
            if tmp is not None:
                container.append(tmp)

def parse_rsids2pos(node):
    res = []
    for i in node['rsids']:
        pos = re.findall('#(\d+)', i)[0]
        res.append(pos)
    return res

def overlap(node1, node2):
    s1 = set(parse_rsids2pos(node1))
    s2 = set(parse_rsids2pos(node2))
    if s1.intersection(s2):
        return True
    else:
        return False
   
    
def compare_node_list(node_list1, node_list2):
    overlap1, overlap2 = [], []

    for i in node_list1:
        tmp = []
        for j in node_list2:
            if overlap(i,j):
                tmp.append(j['name'])
        overlap1.append(tmp)

    for i in node_list2:
        tmp = []
        for j in node_list1:
            if overlap(i,j):
                tmp.append(j['name'])
        overlap2.append(tmp)

    eval_table1 = pd.DataFrame([i['name'] for i in node_list1], columns=['node on manual tree'])

    eval_table1['overlap with auto tree node'] = overlap1
    eval_table2 = pd.DataFrame([i['name'] for i in node_list2], columns=['node on auto tree'])

    eval_table2['overlap with manual tree node'] = overlap2

    eval_table2
    return eval_table1,eval_table2

def compare_tree(tree1, tree2):
    node_list1, node_list2 = [], []
    flatten_tree(node_list1, tree1)
    flatten_tree(node_list2, tree2)
    eval_table1, eval_table2 = compare_node_list(node_list1, node_list2)
    return eval_table1

    

def match_eval(eval_table):
    eval1 = eval_table
    condition = []
    for i in range(eval1.shape[0]):
        item = eval1.iloc[i, :]
        l , r = item['node on manual tree'], item['overlap with auto tree node']
        if len(r) == 1:
            pos1 = pd.DataFrame(parse_rsids(search_node(str(l), manual_tree)), columns=['pos', 'ref', 'alt']).pos.tolist()
            pos2 = pd.DataFrame(parse_rsids(search_node(str(r[0]), tree)), columns=['pos', 'ref', 'alt']).pos.tolist()
            report = generate_report(pos1, pos2)
            if report.shape[0] < 1:
                condition.append('SNP完全匹配')
                continue
        condition.append('不一致')
    eval1['状态'] = condition
    return eval1


def vip_nodes(chipTree, manualTree, frameTree):
    chip_node, manual_node, frame_node = [],[],[]
    flatten_tree(chip_node, chipTree)
    flatten_tree(manual_node, manualTree)
    flatten_tree(frame_node, frameTree)

    eval1, eval2 = compare_node_list(manual_node, chip_node)

    overlap_name = []
    for i in eval1['overlap with auto tree node'].tolist():
        overlap_name+=i

    overlap_nodes = []
    for i in chip_node:
        if i['name'] in overlap_name:
            overlap_nodes.append(i)

    tmpeval1, tmpeval2 = compare_node_list(overlap_nodes, manual_node)

    tmpeval1.columns=['chip_tree_node', 'manual_tree_node']

    eval1, eval2 = compare_node_list(overlap_nodes, frame_node)

    tmp = eval1
    tmp.columns = ['chip_tree_node', 'frame_tree_node']
    chip_tree_node_name = []
    for i in tmp.chip_tree_node.tolist():
        chip_tree_node_name.append(search_node(i, chipTree)['shortName'])
    tmp['chip_tree_node_name'] = chip_tree_node_name
    tmp[['chip_tree_node', 'chip_tree_node_name', 'frame_tree_node']]
    res = pd.merge(tmp, tmpeval1, on='chip_tree_node', how='left')
    res = res[['chip_tree_node', 'chip_tree_node_name', 'manual_tree_node', 'frame_tree_node']]
    return res

def search_pos_in_vip_table(pos, table):
    pos = search_pos(pos, manual_tree)
    f = lambda x: pos in x
    result = table[table['manual_tree_node'].map(f)]
    return result


def out_put_snps(pos, tree):
    tmp = parse_rsids(search_node(pos, tree))
    for i in tmp:
        print(i[0], i[1], i[2])

def compose_node(node_name, tree):
    res = pd.DataFrame(parse_rsids(search_node(node_name, tree)), columns=['pos','ref','alt'])
    res['name'] = node_name
    return res


def zero_match(x):
    return len(x) < 1

def one_match(x):
    return len(x) == 1

def none_zero_match(x):
    return len(x) > 0

def multi_match(x):
    return len(x)>1


class Node(): 
    def __init__(self, node):
        self.node = node
        self.node_pos_set = set(self.node.pos)
    
    def union(self, other):
        new_node = pd.concat((self.node, other.node))
        new_node.drop_duplicates('pos', inplace=True)
        return Node(new_node)
    
    def intersection(self, other):
        new_node = pd.concat((self.node, other.node))
        new_node.drop_duplicates('pos', inplace=True)
        intersection = self.node_pos_set.intersection(other.node_pos_set)
        new_node = new_node[new_node['pos'].isin(intersection)]
        return Node(new_node)
    
    def difference(self, other):
        new_node = pd.concat((self.node, other.node))
        new_node.drop_duplicates('pos', inplace=True)
        difference = self.node_pos_set.difference(other.node_pos_set)
        new_node = new_node[new_node['pos'].isin(difference)]
        return Node(new_node)
    
    def snp_num(self):
        print(self.node.shape[0])
    
    def __repr__(self):
        res = self.node['pos'].astype(str)
        for i in self.node.keys()[1:]:
            res += ',' + self.node[i].astype(str)
#         res = self.node['pos']+' '+self.node['ref']+' '+self.node['alt']
        s = ''
        for i in res.tolist():
            s+=i+'\n'
        return s

def barcode_as_name(tree):
    if 'barcode' in tree.keys():
        tree['name'] = tree['barcode']
    if tree['children'] is None:
        return None
    else:
        for i in tree['children']:
            tmp = barcode_as_name(i)
      
    
def shortName_as_name(tree):
    tree['name'] = tree['shortName']
    if tree['children'] is None:
        return None
    else:
        for i in tree['children']:
            shortName_as_name(i)


def mut_at_node(node_name, barcodes, tree):
    n = Node(compose_node(node_name,tree))
    res = hap_certainty(barcodes, n.node.pos.tolist())

    res['pos'] = res.index

    res.index.name = None

    res = pd.merge(n.node, res, on='pos', how='left')
#     res = res[pd.notna(res[barcode])]
    return res


def get_parent(hapname, tree):
    if tree['children'] is None:
        return None
    for i in tree['children']:
        if i['name'] == hapname:
            return tree['name']
        else:
            tmp = get_parent(hapname, i)
            if tmp is not None:
                return tmp

def chip_to_frame(chipname, chip_tree, frame_tree):
    n = Node(compose_node(chipname, chip_tree))
    match_res = []
    for i in n.node.pos.tolist():
        match_res.append(search_pos(i, frame_tree))
    n.node['frame_node'] = match_res
    return n

def print_row(L):
    s = ''
    for i in L:
        s+=str(i)+','
    s = s[:-1]
    print(s)
    
def print_col(L):
    s = ''
    for i in L:
        s+=str(i)+'\n'
    s = s[:-1]
    print(s)

    
def input_barcodes():
    bars = parse_barcode(input('barcodes:\n'))
    for i in bars:
        if bars.count(i) > 1:
            bars.remove(i)
    return bars


def get_samples(tree):
    jsTree = json.dumps(tree)
    samples = re.findall(bar_pattern, jsTree)
    samples = list(set(samples))
    return samples

def check_sample_location(node_name, bars, shape):
    n = Node(compose_node(node_name, shape.tree))
    for i in shape.get_direct_branches(node_name):
        n_t = Node(compose_node(i, shape.tree))
        n = n.union(n_t)
    res = hap_certainty(bars, n.node.pos.tolist())
    n.node = pd.merge(n.node, res, on='pos', how='left')
    return n


def parse_qc(qc_str):
    qc_str = str(qc_str)
    res = re.findall(qc_pattern, qc_str)
    if len(res) < 1:
        return (0,0)
    res = res[0]
    return eval_code_str(res[0]), eval_code_str(res[1])

def eval_code_str(string):
    if string == 'N':
        return 0
    nums = string.split('-')
    if len(nums) > 1:
        return (int(nums[0])+int(nums[1]))*1.0/2
    else:
        return (int(nums[0].split('+')[0]))

def mut_eval_on_key(n, key):
    n.node.fillna('(N|N)', inplace=True)
    n.node['mut'] = n.node[key].map(parse_qc).map(lambda x:x[0]<x[1])
    return n


def gt_eval(refN, altN):
    if refN > altN:
        return 'ref'
    elif refN < altN:
        return 'alt'
    elif (refN == altN) & (refN>0):
        return 'heterotype'
    else:
        return 'nocall'
    
    
def qc_gt_eval(x):
    tmp = parse_qc(x)
    return gt_eval(tmp[0], tmp[1])    


def compose_ancestry_str(barcode, df):
    sub_df = df[df['barcode'] == barcode]
    if sub_df.shape[0] < 1:
        return '暂无信息'
    ancestry_str = sub_df['province']+ ' '+sub_df['county']
    if pd.isna(ancestry_str.iloc[0]):
        ancestry_str = sub_df['birthplace']
    if pd.isna(sub_df['surname'].iloc[0]):
        return ancestry_str.iloc[0]
    else:
        ancestry_str = ancestry_str+' '+ sub_df['surname']
        return ancestry_str.iloc[0]

def add_ancestry_str(tree, df):
    if 'barcode' in tree.keys():
        tree['ancestry'] = compose_ancestry_str(tree['barcode'], df)
    else:
        for i in tree['children']:
            add_ancestry_str(i)
            
            
def valid_mean(series):
    high = series.mean()+1.5*series.std()
    low = series.mean()-1.5*series.std()
    return series[(series<=high) & (series>=low)].mean()



def compose_ancestry_str(barcode, df):
    sub_df = df[df['barcode'] == barcode]
    if sub_df.shape[0] < 1:
        return '暂无信息'
    if pd.isna(sub_df[['province', 'city', 'county']]).sum(axis=1).iloc[0] == 0 :
        ancestry_str = sub_df['birthplace']
    else:
        sub_df.fillna('', inplace=True)
        ancestry_str = sub_df['province']+' '+ sub_df['city']+ ' ' + sub_df['county']
    if pd.isna(sub_df['surname'].iloc[0]) == False:
        ancestry_str = ancestry_str+' '+ sub_df['surname']
    if pd.isna(sub_df['nation'].iloc[0]) == False:
        ancestry_str =  ancestry_str+' '+sub_df['nation']
    return ancestry_str.iloc[0]

def add_ancestry_str(tree, df):
    if 'barcode' in tree.keys():
        tree['ancestry'] = compose_ancestry_str(tree['barcode'], df)
    else:
        for i in tree['children']:
            add_ancestry_str(i, df)
            
            
def force_submit(handler):
    s = '''插入节点\ttmp\t\t\t\t'''
    s2 = '''删除节点\ttmp\t是\t\t\t'''
    handler.run_excel(s)
    handler.run_excel(s2)
    handler.submit()
    
    
def visualize(pos, bars, width):
    pos = [str(i) for i in range(int(pos)-int(width/2), int(pos)+int(width/2))]
    tmp = hap_certainty([bars], pos)
    qc = tmp[bars].map(parse_qc)
    ref, alt =[],[]
    for i in qc:
        ref.append(i[0])
        alt.append(i[1])
    plt.bar(pos, ref)
    plt.bar(pos, alt, bottom=ref)
    plt.legend(['ref', 'alt'])
    

def check_samples(bars, pos, node_name):

    tmp = hap_certainty(bars, pos)

    tmp = tmp.applymap(qc_gt_eval)

    result = pd.DataFrame(columns=['barcode', 'parent_node', 'ref', 'alt', 'nocall', 'heterotype'])

    for i in tmp.keys():
        item = tmp[i].value_counts()
        item = [i, node_name, item.get('ref', 0), item.get('alt', 0), item.get('nocall', 0), item.get('heterotype', 0)]
        result.loc[result.shape[0]+1] = item

    return result


def eval_qc_table(tb, node_name):
    tmp = tb
    
    tmp = tmp.applymap(qc_gt_eval)

    result = pd.DataFrame(columns=['barcode', 'parent_node', 'ref', 'alt', 'nocall', 'heterotype'])

    for i in tmp.keys():
        item = tmp[i].value_counts()
        item = [i, node_name, item.get('ref', 0), item.get('alt', 0), item.get('nocall', 0), item.get('heterotype', 0)]
        result.loc[result.shape[0]+1] = item

    return result
    


def check_barcode(check_bar, tree):

    node = search_node(check_bar, tree)['name']

    parent = get_parent(node, tree)

    pos = compose_node(parent, tree).pos.tolist()

    bars = []

    for i in search_node(parent,tree)['children']:
        if 'barcode' in i.keys():
            bars.append(i['barcode'])

    result = check_samples(bars, pos, parent)
    return result


def check_barcodes(bars, tree):
    result = pd.DataFrame()
    for i in bars:
        tmp = check_barcode(i, tree)
        result = result.append(tmp)
    return result


def parse_new_sample_pre_hap(s):
    s = s.split(' ')

    s = [i.split('\t') for i in s]

    s = pd.DataFrame(s, columns=['barcode', 'yhaplogroup'])
    s['group'] = s['yhaplogroup'].str[0]
    return s


def name_sample(bar, shape):
    node = shape.search_barcode(bar)
    
    if node is None:
        print('%s not on tree'%bar)
        return ''

    names = shape.parse_pos_full(node)

    names = names[names.pos_name != '']

    if names.shape[0] < 1:
        name = ''
        return name
    
    MF = names[names.pos_name.str[:2] == 'MF']

    if MF.shape[0] > 0:
        name = MF.pos_name.iloc[0]

    else: 
        name = names.pos_name.iloc[0]
    return name


def all_sample_path(tree, res=[], parent_path=['ROOT']):
    name = tree['name']
    if 'barcode' in tree.keys():
        barcode = tree['barcode']
        res.append([barcode, parent_path])
    if 'children' in tree.keys():
        if len(tree['children']) > 0:
            for i in tree['children']:
                get_bar_path(i, res, parent_path+[i['name']])
                
                
def is_chip_user(s):
    return len(re.findall('111-16\d{2}-\d{4}', s)) > 0

def is_out_user(s):
    return len(re.findall('990-\d{4}-\d{4}', s)) > 0

def get_sample_path(tree, parent=['ROOT']):
    name = tree['name']
    if 'barcode' in tree.keys():
        barcode = tree['barcode']
        res.append([barcode, parent])
    else:
        if 'children' in tree.keys():
            if len(tree['children']) > 0:
                for i in tree['children']:
                    get_sample_path(i, parent+[i['name']])

def parse_input_barcode():
    return parse_barcode(input('粘贴输入barcode:\n'))

                
def samtools_query(pos, barcode):
    print('samtools tview -p chrY:%s  ./%s/%s.chrY.bam'%(pos, barcode, barcode))

    
def bam_file_path():
    print('/data6/public/yqaunxudata/YQUANXU/kaijun 111')
    print('/data7/NGS/01.Yfull/01.Ybam/ 990')

    
def get_3x_pos(handler, node_name):
    '''handler, node_name'''
    node, parent = handler.find_node(node_name)
    pos = handler.calc.user_loci_num(node["barcode"], node["rsids"], print_left=True)
    return pos

def generate_qc_table(bars, pos):
    '''
    return (raw_table, eval_table, description_table)
    '''

    tmp1 = hap_certainty(bars, pos)
    try:
        tmp = tmp1.applymap(qc_gt_eval)
    except:
        tmp = None
    
    try:
        result = pd.DataFrame(columns=['barcode', 'parent_node', 'ref', 'alt', 'nocall', 'heterotype'])

        for i in tmp.keys():
            item = tmp[i].value_counts()
            item = [i, tmp, item.get('ref', 0), item.get('alt', 0), item.get('nocall', 0), item.get('heterotype', 0)]
            result.loc[result.shape[0]+1] = item
    except:
        result = None

    return (tmp1, tmp, result)

def formalize_set(_set):
    new_set = [str(i) for i in _set]
    return set(new_set)

def set_comparison(s1, s2):
    s1 = formalize_set(s1)
    s2 = formalize_set(s2)
    diff1 = s1.difference(s2)
    diff2 = s2.difference(s1)
    inner = s1.intersection(s2)
    max_row = max(len(diff1), len(diff2), len(inner))
    result = pd.DataFrame()
    result['%s_common'%(len(inner))] = list(inner) + ([None]*(max_row-len(inner)))
    result['%s_set1_unique'%(len(diff1))] = list(diff1) + ([None]*(max_row-len(diff1)))
    result['%s_set2_unique'%(len(diff2))] = list(diff2) + ([None]*(max_row-len(diff2)))
    return result

def percentage_format(numf):
    return str(round(numf*100, 4))+' %'

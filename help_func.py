import pandas as pd

def taiwan_as_one_province(city_area_info, city_population, all_data, male_data):

    tmp = city_area_info[city_area_info['province'] == '台湾']
    city_area_info= city_area_info[city_area_info['province'] != '台湾']

    merged = ['台湾','台湾省']+tmp[['area', 'population']].sum().tolist()
    city_area_info.loc[city_area_info.index[-1]+1] = merged

    tmp = city_population[city_population['province'] == '台湾']
    city_population= city_population[city_population['province'] != '台湾']

    merged = ['台湾','台湾省']+tmp[['population','weight']].sum().tolist()

    city_population.loc[city_population.index[-1]+1] = merged

    all_data['city'][all_data['province'] == '台湾'] = '台湾省'
    male_data['city'][male_data['province'] == '台湾'] = '台湾省'
    
    return city_area_info, city_population, all_data, male_data
def merge_metrix_info(original_df, patch_df, index_key='barcode'):
    '''
    @param original_df: original pandas dataframe
    @param patch_df: extra information
    @param index_key: unique key to merge original_df and patch_df
    @return: keep information in both original_df and patch_df
    '''
    fix_key_set = set(patch_df.keys())
    original_df_key_set = set(original_df.keys())
    inner_key_set = fix_key_set.intersection(original_df_key_set)
    patch_df = patch_df[list(inner_key_set)]
    fix_part = original_df[original_df[index_key].isin(patch_df[index_key].tolist())]
    original_df = original_df[original_df[index_key].isin(fix_part[index_key])==False]

    patch_df.index = patch_df[index_key]
    fix_part.index = fix_part[index_key]

    fix_mask = patch_df.where(lambda x: pd.isna(x), 1)
    fix_mask = fix_mask.fillna(0)

    data_mask = fix_part.where(lambda x: pd.isna(x), 1)
    data_mask = data_mask.fillna(0)

    mask = (fix_mask + data_mask)[fix_part.keys()]
    mask.where(lambda x: x!=2, None, inplace=True)
    mask.where(lambda x: x>0, None, inplace=True)
    mask.where(lambda x: x!=1, '', inplace=True)

    fix_val = patch_df + mask
    fix_val = fix_val.where(lambda x: pd.notna(x), '')
    fix_part = fix_part.where(lambda x: pd.notna(x), '')
    fixed_val = fix_val + fix_part.astype(str)
    original_df = pd.concat([original_df, fixed_val])
    original_df = original_df[original_df[index_key].notna()]
    return original_df

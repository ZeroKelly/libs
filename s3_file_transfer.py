import os

def upload_to_s3_ancestry(*args):
    '''@param *args: 文件夹名称
    把本地文件夹上传至s3'''
    for name in args:
        cmd = '/home/jovyan/.local/bin/aws s3 cp --recursive ./%s/ s3://genodig-prod/workspaces/users/ancestry/%s/'%(name, name)
        if os.system(cmd) != 0:
            print('fail to excute:\n', cmd)
        else:
            print(cmd, '\n succesful. \n')

def download_from_s3_ancestry(*args):
    '''@param *args: 文件夹名称
    从s3下载文件夹覆盖本地文件夹'''
    for name in args:
        cmd = '/home/jovyan/.local/bin/aws s3 cp --recursive s3://genodig-prod/workspaces/users/ancestry/%s/ ./%s/'%(name, name)
        if os.system(cmd) != 0:
            print('fail to excute:\n', cmd)
        else:
            print(cmd, '\n succesful. \n')

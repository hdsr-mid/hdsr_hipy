# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 08:59:35 2021

@author: hurkmans
"""


api_url = "https://api.dataplatform.knmi.nl/open-data"

if api_key is None: 
    api_key = 'eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6ImNjOWE2YjM3ZjVhODQwMDZiMWIzZGIzZDRjYzVjODFiIiwiaCI6Im11cm11cjEyOCJ9'
if max_keys is None:
    max_keys = 10
if start_after is None:
    start_after = ' '


session = requests.Session()
session.headers.update({"Authorization": api_key})

#start_after_filename_prefix = 'INTER_OPER_R___EV24_H__L3__19110515T000000_19110516T000000_0001'
list_files_response =session.get(f"{api_url}/v1/datasets/{dataset_name}/versions/{dataset_version}/files",
                                   headers={"Authorization": api_key},
                                   params={"maxKeys": max_keys,
                                           "startAfterFilename": start_after})
list_files = list_files_response.json()
dataset_files = list_files.get("files")

temp_path = download_path / 'temp'
temp_path.mkdir(parents=True, exist_ok=True)
    
# Retrieve first file in the list files response
for dataset_file in dataset_files:        
    if not os.path.exists(temp_path / dataset_file['filename']):
        logger.info('Downloading '+dataset_file['filename'])
        filename = dataset_file.get("filename")
        endpoint = f"{api_url}/datasets/{dataset_name}/versions/{dataset_version}/files/{filename}/url"
        #get_file_response = requests.get(endpoint, headers={"Authorization": api_key})
        get_file_response = session.get(endpoint)
        download_url = get_file_response.json().get("temporaryDownloadUrl")
        dataset_file_response = requests.get(download_url)

        # Write dataset file to disk
        p = temp_path / filename
        p.write_bytes(dataset_file_response.content)
    else: 
        logger.info(str(temp_path / dataset_file['filename'])+' exists.')
        
        
        
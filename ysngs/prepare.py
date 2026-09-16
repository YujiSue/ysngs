import os
import json
import datetime
import pandas as pd

def preparInput(script, prop):
    if script == 'scrnaseq':
        return prepare_input_scrnaseq(prop)
    else:
        print(f"Error: Unsupported script '{script}'")
        return None

# scRNA-seq input preparation
def prepare_input_scrnaseq(prop):
    ## Load sample info
    dir = os.path.join(os.environ['HYM_DATA'], prop['dir'])
    df = pd.read_csv(os.path.join(dir, prop['sample_list']), header=0)
    ############################################################
    ## Template:                                              ##
    ##    sample_template(scRNA).csv                          ##
    ############################################################
    ## Set input path
    now = datetime.datetime.now()
    input_data_path = os.path.join(os.environ['HYM_LOG'], f"scrna_{now.strftime('%Y-%m-%d')}.json")
    input_data = []
    ##
    cr_cfg = {}
    for idx,row in df.iterrows():
        if row['SampleID'] not in cr_cfg:
            cr_cfg[row['SampleID']] = []
        cr_cfg[row['SampleID']].append({
            'path' : os.path.join(dir, row['Files']), 'type' : row['Type'] 
        })
    ##
    for sid in cr_cfg:
        if len(cr_cfg[sid]) > 1:
            ## Multi mode
            cfg_path = os.path.join(dir, f'{sid}_config.csv')
            with open(cfg_path, 'w') as f:
                types = [item["type"] for item in cr_cfg[sid]]
                ####
                if 'gex' in types:
                    f.write('[gene-expression]\n')
                    f.write(f"reference,{os.path.join(os.environ['HYM_REF'], prop['gex_reference'])}\n")
                    f.write(f"create-bam,{'true' if prop['export_bam'] else 'false'}\n\n")
                ####
                if 'vdj' in types:
                    f.write('[vdj]\n')
                    f.write(f"reference,{os.path.join(os.environ['HYM_REF'], prop['vdj_reference'])}\n\n")
                ####
                f.write('[libraries]\nfastq_id,fastqs,feature_types\n')
                for sample in cr_cfg[sid]:
                    if sample['type'] == 'gex':
                        f.write(f"{sid},{sample['path']},Gene Expression\n")
                    elif sample['type'].startswith('vdj'):
                            f.write(f"{sid},{sample['path']},upper({sample['type']})\n")
            ##
            input_data.append({
                "sctranscriptome.multi": True,
                "sctranscriptome.config": cfg_path,
                "sctranscriptome.dir": dir,
                "sctranscriptome.out_name": prop['out_name'],
                "sctranscriptome.thread": prop['thread'],
                "sctranscriptome.ram": prop['ram']
                })
        else:
            ## Single mode
            input_data.append({
                "sctranscriptome.multi": False,
                "sctranscriptome.export_bam": prop['export_bam'],
                "sctranscriptome.reference": os.path.join(os.environ['HYM_REF'], prop['gex_reference']),
                "sctranscriptome.dir": dir,
                "sctranscriptome.subdir": cr_cfg[sid][0]['path'],
                "sctranscriptome.file_prefix": sid,
                "sctranscriptome.out_name": prop['out_name'],
                "sctranscriptome.thread": prop['thread'],
                "sctranscriptome.ram": prop['ram']
            })
    ## Save input data to JSON
    json.dump(input_data, open(input_data_path, 'w'))

    ## Return input path
    return input_data_path
        

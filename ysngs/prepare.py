import os
import json
import datetime
import pandas as pd

def makeInputPath(prefix):
  now = datetime.datetime.now()
  return os.path.join(os.environ['HYM_LOG'], f"{prefix}_{now.strftime('%Y-%m-%d')}.json")

def preparInput(script, prop):
  if script == 'dlfq':
    return prepare_input_dlfq(prop)
  if script == 'mkidx':
    return prepare_input_mkidx(prop)
  if script == 'scrnaseq':
    return prepare_input_scrnaseq(prop)
  else:
    print(f"Error: Unsupported script '{script}'")
    return None

# Download data from sequence archive
def prepare_input_dlfq(prop):
    ## Set path
    input_data_path = makeInputPath('dlfq')
    input_data = []
    ## Set IDs
    if 'ids' not in prop:
        prop['ids'] = []
    if 'list' in prop and os.path.exists(prop['list']):
        df = pd.read_csv(os.path.join(os.environ['HYM_DATA'], list), header=None)
        prop['ids'].extend(df.iloc[:, 0].tolist())
    ## Set input 
    for file_id in prop['ids']:
        input_data.append({
            "dlfq.data_id": file_id,
            "dlfq.split": prop['split'],
            "dlfq.out_dir": os.path.join(os.environ['HYM_DATA'], prop['out_dir']),
            "dlfq.thread": prop['thread']
        })
    ## Save input data to JSON
    json.dump(input_data, open(input_data_path, 'w'))
    ## Return input path
    return input_data_path

# Make index file(s)
def prepare_input_mkidx(prop):
    ## Set path
    input_data_path = makeInputPath('dlfq')
    input_data = []
    ## Select apps
    input_data.append({
        "mkidx.use_hts": prop['use_hts'],
        "mkidx.use_bwa": prop['use_bwa'],
        "mkidx.use_bowtie": prop['use_bowtie'],
        "mkidx.use_gatk": prop['use_gatk'],
        "mkidx.use_star": prop['use_star'],
        "mkidx.use_hisat": prop['use_hisat'],
        "mkidx.use_rsem": prop['use_rsem'],
        "mkidx.use_cr": prop['use_cr'],
        "mkidx.out_dir": prop["out_dir"],
        "mkidx.ref_fasta": prop["fasta"],
        "mkidx.ref_gtf": prop["gtf"],
        "mkidx.ref_label": prop["name"],
        "mkidx.mapper_path": prop["app_path"],
        "mkidx.thread": prop["thread"]
    })
    ## Save input data to JSON
    json.dump(input_data, open(input_data_path, 'w'))
    ## Return input path
    return input_data_path

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
        

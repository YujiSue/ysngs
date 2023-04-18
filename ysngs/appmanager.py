from ysngs import installer

APP_LIST = {
  'igv': { 'name': 'IGV', 'ver': '2.13.2', 'checker': installer.checkIGV, 'verchecker': installer.checkVerIGV, 'installer': installer.installIGV },
  'sra': { 'name': 'SRA-toolkit', 'ver': '2.11.3', 'checker': installer.checkSRA, 'verchecker': installer.checkVerSRA, 'installer': installer.installSRA },
  #'root': { 'name': 'CERN ROOT', 'ver': '6.28.02', 'checker': installer.checkROOT, 'verchecker': installer.checkVerROOT, 'installer': installer.installROOT },
  'blast': {},
  'samtools': { 'name': 'SAMTools', 'ver': '1.16', 'checker': installer.checkSamtools, 'verchecker': installer.checkVerSamtools, 'installer': installer.installSamtools },
  'bcftools': { 'name': 'BCFTools', 'ver': '1.16', 'checker': installer.checkBCFtools, 'verchecker': installer.checkVerBCFtools, 'installer': installer.installBCFtools },
  'htslib': { 'name': 'HTSlib', 'ver': '1.16', 'checker': installer.checkHTSlib, 'verchecker': installer.checkVerHTSlib, 'installer': installer.installHTSlib },
  'picard': { 'name': 'Picard', 'ver': '2.27.4', 'checker': installer.checkPicard, 'verchecker': installer.checkVerPicard, 'installer': installer.installPicard },
  'gatk': { 'name': 'GATK', 'ver': '4.2.4.0', 'checker': installer.checkGATK, 'verchecker': installer.checkVerGATK, 'installer': installer.installGATK },
  'cutadapt':{ 'name': 'CutAdapt', 'ver': 'latest', 'checker': installer.checkCut, 'verchecker': installer.checkVerCut, 'installer': installer.installCut },
  'fastqc': { 'name': 'FastQC', 'ver' : '0.11.9', 'checker': installer.checkFQC, 'verchecker': installer.checkVerFQC, 'installer': installer.installFQC },
  'fastp': { 'name': 'fastp', 'ver' : '0.23.2', 'checker': installer.checkFP, 'verchecker': installer.checkVerFP, 'installer': installer.installFP },
  'bwa': { 'name': 'BWA', 'ver': '0.7.17', 'checker': installer.checkBWA, 'verchecker': installer.checkVerBWA, 'installer': installer.installBWA },
  'bowtie2': { 'name': 'Bowtie2', 'ver': '2.4.4', 'checker': installer.checkBowtie, 'verchecker': installer.checkVerBowtie, 'installer': installer.installBowtie },
  'star': { 'name': 'STAR', 'ver': '2.7.9a', 'checker': installer.checkSTAR, 'verchecker': installer.checkVerSTAR, 'installer': installer.installSTAR },
  'hisat2': { 'name': 'HISAT2', 'ver': '2.2.1', 'checker': installer.checkHISAT2, 'verchecker': installer.checkVerHISAT2, 'installer': installer.installHISAT2 },
  'ivc': { 'name': 'Strelka', 'ver': '2.9.2', 'checker': installer.checkStrelka, 'verchecker': installer.checkVerStrelka, 'installer': installer.installStrelka },
  'manta': { 'name': 'Manta', 'ver': '1.6.0', 'checker': installer.checkManta, 'verchecker': installer.checkVerManta, 'installer': installer.installManta },
  'tvc': { 'name': 'TorrentVariantCaller', 'ver': '5.12.1', 'checker': installer.checkTVC, 'verchecker': installer.checkVerTVC, 'installer': installer.installTVC },
  'gdv': { 'name': 'GoogleDeepVariant', 'ver': '1.3.0', 'checker': installer.checkGDV, 'verchecker': installer.checkVerGDV, 'installer': installer.installGDV },
  'cnvnator': { 'name': 'CNVnator', 'ver': 'latest', 'checker': installer.checkCNVnator, 'verchecker': installer.checkVerCNVnator, 'installer': installer.installCNVnator },
  'delly': {'name': 'DELLY', 'ver': 'latest', 'checker': installer.checkDELLY, 'verchecker': installer.checkVerDELLY, 'installer': installer.installDELLY },
  'lumpy': { 'name': 'LUMPY', 'ver': 'latest', 'checker': installer.checkLUMPY, 'verchecker': installer.checkVerLUMPY, 'installer': installer.installLUMPY },
  'evep': { 'name': 'ensembl-vep', 'ver': 'latest', 'checker': installer.checkEVEP, 'verchecker': installer.checkVerEVEP, 'installer': installer.installEVEP },
  'sift': { 'name': 'SIFT-4G', 'ver': 'latest', 'checker': installer.checkSift, 'verchecker': installer.checkVerSift, 'installer': installer.installSift, 'dbinstaller': installer.installSiftDB },
  'interpro': {'name': 'InterProScan', 'ver': '5.60-92.0', 'checker': installer.checkInterPro, 'verchecker': installer.checkVerInterPro, 'installer': installer.installInterPro },
  'hlahd': { 'name': 'HLA-HD', 'ver': '1.5.0', 'checker': installer.checkHLAHD, 'verchecker': installer.checkVerHLAHD, 'installer': installer.installHLAHD },
  'htseq': { 'name': 'HTSeq', 'ver': 'latest', 'checker': installer.checkHTSeq, 'verchecker': installer.checkVerHTSeq, 'installer': installer.installHTSeq },
  'fcount': { 'name': 'featureCounts', 'ver':'latest', 'checker': installer.checkFeatCounts, 'verchecker': installer.checkVerFeatCounts, 'installer': installer.installFeatCounts },
  'cufflink': { 'name': 'Cufflinks', 'ver': '2.2.1', 'checker': installer.checkCuff, 'verchecker': installer.checkVerCuff, 'installer': installer.installCuff },
  'rsem': { 'name': 'RSEM', 'ver': '1.3.3', 'checker': installer.checkRSEM, 'verchecker': installer.checkVerRSEM, 'installer': installer.installRSEM },
  'r-bioc': { 'name': 'BiocManager', 'ver': 'latest', 'checker': installer.checkBM, 'installer': installer.installBM },
  
  'macs2': { 'name': 'MACS2', 'ver':'v2.2.7.1', 'checker': installer.checkMACS, 'verchecker': installer.checkVerMACS, 'installer': installer.installMACS },
  
  'meme': { 'name': 'MEME', 'ver':'5.4.1', 'checker': installer.checkMEME, 'verchecker': installer.checkVerMEME, 'installer': installer.installMEME }
}

class AppManager:
    def __init__(self, config):
        self.cfg = config
    def installable(self):
        info = {}
        keys = APP_LIST.keys()
        for k in keys:
            info[k] = {
                'name': APP_LIST[k]['name'],
                'version': APP_LIST[k]['ver']
            }
        return info
    def install(self, name):
        if name in APP_LIST:
            app_check = APP_LIST[name]['checker'](self.cfg)
            print('Check '+APP_LIST[name]['name']+' ...', 'Installed. (ver. '+APP_LIST[name]['verchecker'](self.cfg)+')' if app_check else 'Not installed.')
            if not app_check:
                APP_LIST[name]['installer'](self.cfg, APP_LIST[name]['ver'])
        else:
            print(name, 'is not installable by this library.')

    # def update(self):
    # def uninstall(self):
    def version(self, name):
        return APP_LIST[name]['verchecker'](self.cfg)
    # def test(name):
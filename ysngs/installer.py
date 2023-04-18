import os
import subprocess
from subprocess import PIPE 
from ysngs import common

# Moirei
def checkMoirei(cfg):
  ret = common.execCmd('which moirei')
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerMoirei(cfg):
  ret = common.execCmd('moirei --version')
  assert ret[0], 'Moirei is not installed.'
  print('> ', ret[1].strip())
def installMoirei(cfg, ver):
  print('Install Moirei ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('git clone ')
  assert common.execCmd('git clone ')
  os.chdir(os.path.join(cfg.TEMPORAL, 'Lemegeton'))
  assert common.execCmd('')
  assert common.execCmd('./Goetia -soai')
  os.chdir(os.path.join(cfg.TEMPORAL, 'Moirei'))
  assert common.execCmd('make -j' + str(cfg.thread), verbose=True)[0], 'Make error.'
  assert common.execCmd('sudo make install', verbose=True)[0], 'Install error.'
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm -r *')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  checkVerMoirei(cfg)

# IGV
def checkIGV(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'igv/igv.sh')) 
def checkVerIGV(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'igv/igv.sh --version'), showcmd=False)
  assert ret[0], 'IGV is not installed.'
  ver = ret[1].splitlines()[-1]
  return ver
def installIGV(cfg, ver):
  print('Install IGV...')
  os.chdir(cfg.TEMPORAL)
  vnum = ver.split('.')
  short_ver = '.'.join([vnum[0], vnum[1]])
  assert common.execCmd('wget https://data.broadinstitute.org/igv/projects/downloads/' + short_ver + '/IGV_Linux_' + ver + '_WithJava.zip')[0], 'Download error.'
  assert common.execCmd('unzip ./IGV_Linux_' + ver + '_WithJava.zip')[0], 'Expansion error.'
  os.makedirs(os.path.join(cfg.APPS_DIR, 'igv'), exist_ok=True)
  assert common.execCmd('mv IGV_Linux_' + ver + '/* ' + os.path.join(cfg.APPS_DIR, 'igv'))[0], 'File move error.'
  assert common.execCmd('rm -r ./IGV*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerIGV(cfg))

# SRA-Toolkit
def checkSRA(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'sra/test-sra')) 
def checkVerSRA(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'sra/test-sra | grep "NCBI SRA Toolkit release version:" | head -n 1'), showcmd=False)
  assert ret[0], 'SRA-Toolkit is not installed.'
  ver = ret[1].replace('NCBI SRA Toolkit release version:','').replace('.<br/>','').strip()
  return ver
def installSRA(cfg, ver):
  print('Install SRA-Toolkit...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/' + ver + '/sratoolkit.'  + ver + '-ubuntu64.tar.gz')[0]
  assert common.execCmd('tar xvf ./sratoolkit.' + ver + '-ubuntu64.tar.gz')[0]
  os.makedirs(os.path.join(cfg.APPS_DIR, 'sra'), exist_ok=True)
  assert common.execCmd('cp -r ./sratoolkit*/bin/* ' + os.path.join(cfg.APPS_DIR, 'sra'))[0]
  assert common.execCmd('rm -r ./sratoolkit*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('> ver.', checkVerSRA(cfg))
  common.execCmd(os.path.join(cfg.APPS_DIR, 'sra', 'vdb-config') + '--interactive')
  
# Samtools
def checkSamtools(cfg):
  ret = common.execCmd('which samtools', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerSamtools(cfg):
  ret = common.execCmd('samtools --version | head -n 1', showcmd=False)
  assert ret[0], 'Samtools is not installed.'
  ver = ret[1].replace('samtools','').strip()
  return ver
def installSamtools(cfg, ver):
  print('Install Samtools...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('curl -LO https://github.com/samtools/samtools/releases/download/' + ver + '/samtools-' + ver + '.tar.bz2')[0], 'Code download error.'
  assert common.execCmd('tar xvf ./samtools-' + ver + '.tar.bz2')[0], 'Expansion error.'
  os.chdir(os.path.join(cfg.TEMPORAL, 'samtools-' + ver))
  assert common.execCmd('./configure', verbose=True)[0], 'Configure error.'
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
  assert common.execCmd('sudo make install', verbose=True)[0], 'Install error.'
  assert common.execCmd('cp '+cfg.TEMPORAL+'/samtools-' + ver + '/examples/ex1.fa '+cfg.REFERENCE_DIR)[0]
  assert common.execCmd('gunzip '+cfg.TEMPORAL+'/samtools-' + ver + '/examples/ex1.sam.gz')[0]
  assert common.execCmd('cp '+cfg.TEMPORAL+'/samtools-' + ver + '/examples/ex1.sam ' +cfg.TEST_DIR)[0]
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm -r ./samtools*')[0]
  os.chdir(cfg.WORK_SPACE)
  assert common.execCmd('samtools view -b -T '+cfg.REFERENCE_DIR+'/ex1.fa -o '+cfg.TEST_DIR+'/ex1.bam '+cfg.TEST_DIR+'/ex1.sam')[0]
  assert common.execCmd('samtools fastq '+cfg.TEST_DIR+'/ex1.bam > '+cfg.TEST_DIR+'/ex1.fq')[0]
  print('Completed.')
  print('> ver.', checkVerSamtools(cfg))

# BCFtools
def checkBCFtools(cfg):
  ret = common.execCmd('which bcftools', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerBCFtools(cfg):
  ret = common.execCmd('bcftools --version | head -n 1', showcmd=False)
  assert ret[0], 'BCFtools is not installed.'
  ver = ret[1].replace('bcftools','').strip()
  return ver
def installBCFtools(cfg, ver):
    print('Install BCFtools ...')
    os.chdir(cfg.TEMPORAL)
    assert common.execCmd('wget https://github.com/samtools/bcftools/releases/download/' + ver + '/bcftools-' + ver + '.tar.bz2')[0], 'Download error.'
    assert common.execCmd('tar xvf ./bcftools-' + ver + '.tar.bz2')[0], 'Expansion error.'
    os.chdir(cfg.TEMPORAL+'/bcftools-' + ver)
    assert common.execCmd('./configure', verbose=True)[0], 'Configure error.'
    assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
    assert common.execCmd('sudo make install', verbose=True)[0], 'Install error.'
    os.chdir(cfg.TEMPORAL)
    assert common.execCmd('rm -r ./bcftools*')[0]
    os.chdir(cfg.WORK_SPACE)
    print('Completed.')
    print('> ver.', checkVerBCFtools(cfg))

# HTSlibrary
def checkHTSlib(cfg):
  return os.path.exists('')
def checkVerHTSlib(cfg):
  return ''
def installHTSlib(cfg, ver):
  print('Install BCFtools ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://github.com/samtools/htslib/releases/download/' + ver + '/htslib-' + ver + '.tar.bz2')[0], 'Download error.'
  assert common.execCmd('tar xvf ./htslib-' + ver + '.tar.bz2')[0], 'Expansion error.'
  os.chdir(cfg.TEMPORAL+'/htslib-' + ver)
  assert common.execCmd('./configure', verbose=True)[0], 'Configure error.'
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
  assert common.execCmd('sudo make install', verbose=True)[0], 'Install error.'
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm -r ./htslib*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  checkVerHTSlib(cfg)

# Picard
def checkPicard(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'picard.jar'))
def checkVerPicard(cfg):
  ret = common.execCmd('java -jar ' + os.path.join(cfg.APPS_DIR, 'picard.jar MarkDuplicates --version'), showcmd=False)
  ver = ret[2].replace('Version:','').strip()
  return ver
def installPicard(cfg, ver):
  print('Install Picard ...')
  assert common.execCmd('wget https://github.com/broadinstitute/picard/releases/download/' + ver + '/picard.jar -O ' + cfg.APPS_DIR + '/picard.jar')[0]
  print('Completed.')
  os.chdir(cfg.WORK_SPACE)
  print('>ver.', checkVerPicard(cfg))
  
# GATK
def checkGATK(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'gatk'))
def checkVerGATK(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'gatk/gatk --list'), showcmd=False)
  ver = ret[2].splitlines()[0].split('-')[-2]
  return ver
def installGATK(cfg, ver):
  print('Install GATK ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://github.com/broadinstitute/gatk/releases/download/' + ver +'/gatk-' + ver + '.zip')[0], 'Download error.'
  assert common.execCmd('unzip -o ./gatk-' + ver + '.zip')[0], 'Expansion error.'
  os.makedirs(os.path.join(cfg.APPS_DIR, 'gatk'), exist_ok=True)
  assert common.execCmd('mv ./gatk-' + ver + '/* ' + os.path.join(cfg.APPS_DIR, 'gatk'))[0]
  assert common.execCmd('rm -r '+cfg.TEMPORAL+'/gatk*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerGATK(cfg))
  
# CutAdapt
def checkCut(cfg):
  ret = common.execCmd('which cutadapt', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerCut(cfg):
  ret = common.execCmd('cutadapt --version', showcmd=False)
  ver = ret[1].strip()
  return ver
def installCut(cfg, ver):
  print('Install cutadapt ...')
  assert common.execCmd('pip install cutadapt', verbose=True)[0], 'Install error.'
  print('Completed.')
  print('>ver.', checkVerCut(cfg))

# FastQC
def checkFQC(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'FastQC/fastqc'))
def checkVerFQC(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'FastQC/fastqc --version &'), showcmd=False)
  assert ret[0], 'fastqc is not installed.'
  ver = ret[1].replace('FastQC','').strip()
  return ver
def installFQC(cfg, ver):
  print('Install FastQC  ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v' + ver + '.zip')[0], 'Download error.'
  assert common.execCmd('unzip fastqc_v' + ver + '.zip')[0], 'Expansion error.'
  assert common.execCmd('mv FastQC ' + cfg.APPS_DIR)[0], 'File move error.'
  assert common.execCmd('chmod a+x ' + os.path.join(cfg.APPS_DIR, 'FastQC/fastqc'))[0], 'Permission setting error.'
  assert common.execCmd('rm fastqc_v' + ver + '.zip')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerFQC(cfg))
  
  
# fastp
def checkFP(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'fastp'))
def checkVerFP(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'fastp --version'), showcmd=False) 
  ver = ret[2].replace('fastp','').strip()
  return ver
def installFP(cfg, ver):
  print('Install fastp  ...')
  assert common.execCmd('wget http://opengene.org/fastp/fastp -O ' + os.path.join(cfg.APPS_DIR, 'fastp'))[0], 'Download error.'
  assert common.execCmd('chmod a+x ' + os.path.join(cfg.APPS_DIR, 'fastp'))[0], 'Permission setting error.'
  print('Completed.')
  print('>ver.', checkVerFP(cfg))

# BWA
def checkBWA(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'bwa'))
def checkVerBWA(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'bwa 2>&1 >/dev/null | grep Version:'), showcmd=False)
  ver = ret[1].replace('Version:', '').strip()
  return ver
def installBWA(cfg, ver):
  print('Install BWA  ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://jaist.dl.sourceforge.net/project/bio-bwa/bwa-' + ver+'.tar.bz2')[0], 'Download error.'
  assert common.execCmd('tar xvf ./bwa-' + ver +'.tar.bz2')[0], 'Expansion error.'
  os.chdir(os.path.join(cfg.TEMPORAL, 'bwa-' + ver))
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
  assert common.execCmd('cp bwa ' + cfg.APPS_DIR)[0], 'File copy error.'
  assert common.execCmd('rm -r ' + os.path.join(cfg.TEMPORAL, 'bwa*'))[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerBWA(cfg))
  
# Bowtie2
def checkBowtie(cfg):
  ret = common.execCmd('which bowtie2', showcmd=False)
  return os.path.exists(ret[1].strip())
def checkVerBowtie(cfg):
  ret = common.execCmd('bowtie2 --version | head -n 1', showcmd=False)
  assert ret[0], 'Bowtie2 is not installed.'
  ver = ret[1].split()[-1].strip()
  return ver
def installBowtie(cfg, ver):
  print('Install Bowtie2 ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://jaist.dl.sourceforge.net/project/bowtie-bio/bowtie2/' + ver + '/bowtie2-' + ver + '-source.zip')[0], 'Download error.'
  assert common.execCmd('unzip -o ./bowtie2-' + ver + '-source.zip')[0], 'Expansion  error.'
  os.chdir(os.path.join(cfg.TEMPORAL, 'bowtie2-' + ver))
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
  assert common.execCmd('sudo make install', verbose=True)[0], 'Install error.'
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm -r ./bowtie2*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerBowtie(cfg))

# STAR
def checkSTAR(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'STAR'))
def checkVerSTAR(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'STAR | grep version'), showcmd=False)
  ver = ret[1].split('=')[-1].strip()
  return ver
def installSTAR(cfg, ver):
  print('Install STAR ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://github.com/alexdobin/STAR/archive/' + ver + '.tar.gz')[0], 'Download error.'
  assert common.execCmd('tar xvf ./' + ver + '.tar.gz')[0], 'Expansion error.'
  os.chdir(os.path.join(cfg.TEMPORAL, 'STAR-' + ver + '/source'))
  assert common.execCmd('make -j8 STAR', verbose=True)[0], 'Make error.'
  assert common.execCmd('cp STAR ' + cfg.APPS_DIR)[0]
  assert common.execCmd('rm -r ' + os.path.join(cfg.TEMPORAL, 'STAR*'))[0]
  assert common.execCmd('rm ' + os.path.join(cfg.TEMPORAL, ver + '*'))[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerSTAR(cfg))
  
# HISAT2
def checkHISAT2(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'hisat2', 'hisat2'))
def checkVerHISAT2(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'hisat2', 'hisat2 --version | head -n 1'), showcmd=False)
  ver = ret[1].split(' ')[-1].strip()
  return ver
def installHISAT2(cfg, ver):
  print('Install HISAT2...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('curl -o hisat2.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download')[0], 'Download error.'
  assert common.execCmd('unzip hisat2.zip')[0], 'Expansion error.'
  assert common.execCmd('mv hisat2-' + ver + ' ' + os.path.join(cfg.APPS_DIR, 'hisat2'))[0]
  assert common.execCmd('rm hisat2.zip')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerHISAT2(cfg))

# IVC
def checkStrelka(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'strelka'))
def checkVerStrelka(cfg):
  ret = common.execCmd('python2 '+os.path.join(cfg.APPS_DIR, 'strelka/bin/configureStrelkaGermlineWorkflow.py --version'), showcmd=False)
  assert ret[0], 'Strelka is not installed.'
  ver = ret[1].strip()
  return ver
#def installIsaac(cfg, ver):
def installStrelka(cfg, ver):
  print('Install Strelka...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://github.com/Illumina/strelka/releases/download/v' + ver + '/strelka-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Download error.'
  assert common.execCmd('tar xvjf strelka-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Expansion error.'
  assert common.execCmd('mv -f strelka-' + ver + '.centos6_x86_64' + ' ' + os.path.join(cfg.APPS_DIR, 'strelka'))[0]
  assert common.execCmd('rm -r strelka*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerStrelka(cfg))
 
# Manta
def checkManta(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'Manta/bin/configManta.py'))
def checkVerManta(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'Manta/bin/configManta.py --version'), showcmd=False)
  assert ret[0], 'Manta is not installed.'
  return ret[1].strip()
def installManta(cfg, ver):
  print('Install Manta...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://github.com/Illumina/manta/releases/download/v' + ver + '/manta-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Download error.'
  assert common.execCmd('tar xvjf manta-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Expansion error.'
  assert common.execCmd('mv ./manta-' + ver + '.centos6_x86_64 ' + os.path.join(cfg.APPS_DIR, 'Manta'))[0]
  assert common.execCmd('rm -r ./manta*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerManta(cfg))

# TVC
def checkTVC(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'TVC/bin/tvc'))
def checkVerTVC(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'TVC/bin/tvc') + ' --version')
  assert ret[0], 'Torrent VariantCaller is not installed.'
  ver = ret[1].split(' ')[1]
  return ver
def installArmadillo4TVC(dir):
  os.chdir(dir)
  assert common.execCmd('wget http://updates.iontorrent.com/updates/software/external/armadillo-4.600.1.tar.gz')[0], 'Download error.'
  assert common.execCmd('tar xvzf armadillo-4.600.1.tar.gz')[0]
  os.chdir('armadillo-4.600.1')
  assert common.execCmd("sed -i 's:^// #define ARMA_USE_LAPACK$:#define ARMA_USE_LAPACK:g' include/armadillo_bits/config.hpp")[0]
  assert common.execCmd("sed -i 's:^// #define ARMA_USE_BLAS$:#define ARMA_USE_BLAS:g'   include/armadillo_bits/config.hpp")[0]
  assert common.execCmd('cmake .', verbose=True)[0], 'CMake error.'
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
def installBamtools4TVC(dir):
  os.chdir(dir)
  assert common.execCmd('wget updates.iontorrent.com/updates/software/external/bamtools-2.4.0.20150702+git15eadb925f.tar.gz')[0], 'Download error.'
  assert common.execCmd('tar xvzf bamtools-2.4.0.20150702+git15eadb925f.tar.gz')[0], 'Expansion error.'
  os.makedirs('bamtools-2.4.0.20150702+git15eadb925f-build',exist_ok=True)
  os.chdir('bamtools-2.4.0.20150702+git15eadb925f-build')
  assert common.execCmd('cmake ../bamtools-2.4.0.20150702+git15eadb925f -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo', verbose=True)[0], 'CMake error.'
  common.execCmd('make -j8', verbose=True)
def installHTSlib4TVC(dir):
  os.chdir(dir)
  assert common.execCmd('wget --no-check-certificate https://github.com/samtools/htslib/archive/1.2.1.tar.gz -O htslib-1.2.1.tar.gz')[0], 'Download error.'
  assert common.execCmd('tar xvzf htslib-1.2.1.tar.gz')[0], 'Expansion error.'
  assert common.execCmd('ln -s htslib-1.2.1 htslib')[0]
  os.chdir('htslib-1.2.1')
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
def installSamtools4TVC(dir, dest):
  os.chdir(dir)
  assert common.execCmd('wget --no-check-certificate https://github.com/samtools/samtools/archive/1.2.tar.gz -O samtools-1.2.tar.gz')[0], 'Download error.'
  assert common.execCmd('tar xvzf samtools-1.2.tar.gz')[0], 'Expansion error.'
  os.chdir('samtools-1.2')
  assert common.execCmd('make -j8', verbose=True)[0]
  assert common.execCmd('cp ./samtools ' + dest)[0]
def installTVC(cfg, ver):
  print('Install Torrent VariantCaller...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget http://updates.iontorrent.com/tvc_standalone/tvc-' + ver+'.tar.gz')[0], 'Download error.'
  build = os.path.join(cfg.TEMPORAL, 'temp')
  os.makedirs(build, exist_ok=True)
  assert common.execCmd('cp '+cfg.TEMPORAL+'/tvc-' + ver + '.tar.gz ' + build + '/tvc-' + ver + '.tar.gz')[0]
  install = os.path.join(cfg.APPS_DIR, 'TVC')
  os.makedirs(install, exist_ok=True)
  os.makedirs(install + '/bin', exist_ok=True)
  installArmadillo4TVC(build)
  installBamtools4TVC(build)
  installHTSlib4TVC(build)
  installSamtools4TVC(build, install + '/bin/')
  os.chdir(build)
  assert common.execCmd('tar xvzf ./tvc-' + ver + '.tar.gz')[0]
  os.makedirs(build+'/build',exist_ok=True)
  os.chdir(build+'/build')
  assert common.execCmd('cmake '+build+'/tvc-' + ver + ' -DCMAKE_INSTALL_PREFIX:PATH=' + install + ' -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo', verbose=True)[0], 'CMake error.'
  assert common.execCmd('make -j8 install', verbose=True)[0], 'Install error.'
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/examples/example1/reference* ' + cfg.REFERENCE_DIR)[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/examples/example1/test* '+cfg.TEST_DIR)[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/pluginMedia/configs/description.json '+cfg.PREFERENCE_DIR)[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/pluginMedia/parameter_sets/ampliseq_germline_lowstringency_pgm_parameters.json ' + os.path.join(cfg.PREFERENCE_DIR, 'germline_low_stringency.json'))[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/pluginMedia/parameter_sets/ampliseq_germline_lowstringency_p1_parameters.json ' + os.path.join(cfg.PREFERENCE_DIR, 'germline_low_stringency_proton.json'))[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/pluginMedia/parameter_sets/targetseq_germline_lowstringency_p1_parameters.json ' + os.path.join(cfg.PREFERENCE_DIR, 'germline_low_stringency_targetseq.json'))[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/pluginMedia/parameter_sets/ampliseq_somatic_lowstringency_pgm_parameters.json ' + os.path.join(cfg.PREFERENCE_DIR, 'somatic_low_stringency.json'))[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/pluginMedia/parameter_sets/ampliseq_somatic_lowstringency_p1_parameters.json ' + os.path.join(cfg.PREFERENCE_DIR, 'somatic_low_stringency_proton.json'))[0]
  assert common.execCmd('cp '+build+'/tvc-' + ver + '/share/TVC/sse/* '+cfg.PREFERENCE_DIR)[0]
  assert common.execCmd(cfg.APPS_DIR+'/TVC/bin/samtools bam2fq '+cfg.TEST_DIR+'/test.bam > '+cfg.TEST_DIR+'/test.fq')[0]
  assert common.execCmd('rm -r '+build)[0]
  assert common.execCmd('rm '+cfg.TEMPORAL+'/tvc*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerTVC(cfg))

# GDV
def checkGDV(cfg):
  ret = common.execCmd('docker images')
  return ret[0] and 'google/deepvariant' in ret[1]
def checkVerGDV(cfg):
  ret = common.execCmd('docker images | grep google/deepvariant')
  assert ret[0], 'Google DeepVariant is not installed.'
  ver = ret[1].split()[1]
  return ver
def installGDV(cfg, ver):
  print('Install Google DeepVariant...')
  assert common.execCmd('sudo docker pull google/deepvariant:'+ver)[0]
  print('Completed.')
  print('>ver.', checkVerGDV(cfg))

# CNVnator
def checkCNVnator(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'CNVnator', 'cnvnator'))
def checkVerCNVnator(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'CNVnator', 'cnvnator 2>&1 >/dev/null | grep CNVnator'), showcmd=False)
  assert ret[0], 'CNVnator is not installed.'
  ver = ret[1].split()[1].strip()
  return ver
def installCNVnator(cfg, ver):
  print('Install CNVnator...')
  assert common.execCmd('echo deb http://security.ubuntu.com/ubuntu xenial-security main | sudo tee -a /etc/apt/sources.list')[0]
  assert common.execCmd('sudo apt-get update')[0]
  assert common.execCmd('sudo apt-get install libssl1.0.0 libreadline-dev')[0]
  os.chdir(cfg.APPS_DIR)
  assert common.execCmd('git clone --recursive https://github.com/abyzovlab/CNVnator.git')[0], 'Download error.'
  os.chdir(os.path.join(cfg.APPS_DIR, 'CNVnator'))
  assert common.download('https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2')[0]
  assert common.execCmd('tar -xvf samtools-1.2.tar.bz2')[0]
  assert common.execCmd('rm samtools-1.2.tar.bz2')[0]
  os.chdir(os.path.join(cfg.APPS_DIR, 'CNVnator', 'samtools-1.2'))
  assert common.execCmd('make -j8')[0], 'Samtools make failed.'
  os.chdir(os.path.join(cfg.APPS_DIR, 'CNVnator'))
  assert common.execCmd('ln -s samtools-1.2 samtools')[0]
  os.environ['ROOTSYS'] = os.path.join(cfg.APPS_DIR, 'root')
  assert common.execCmd('make -j8')[0], 'Make failed.'
  os.chdir(os.path.join(cfg.WORK_SPACE))
  print('Completed.')
  print('>ver.', checkVerCNVnator(cfg))

# DELLY
def checkDELLY(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'delly/src/delly'))
def checkVerDELLY(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR ,'delly/src/delly 2>&1 >/dev/null | grep Version:'), showcmd=False)
  ver = ret[1].split()[-1][:-1]
  return ver
def installDELLY(cfg, ver):
  print('Install DELLY...')
  os.chdir(cfg.APPS_DIR)
  assert common.execCmd('git clone --recursive https://github.com/dellytools/delly.git')[0], 'Download error.'
  os.chdir('./delly')
  assert common.execCmd('sudo make all -j8')[0], 'Compile error.'
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerDELLY(cfg))


# LUMPY
def checkLUMPY(cfg):
  ret = common.execCmd('which lumpy', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerLUMPY(cfg):
  ret = common.execCmd('lumpy 2>&1 >/dev/null | grep Program', showcmd=False)
  ver = ret[1].split()[-1][1:-1]
  return ver
def installLUMPY(cfg, ver):
  print('Install LUMPY...')
  ret = common.execCmd('which autoreconf')
  if not (ret[0] and os.path.exists(ret[1].strip())):
    assert common.execCmd('sudo apt-get install automake')[0]
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('git clone --recursive https://github.com/arq5x/lumpy-sv.git')[0], 'Download error.'
  os.chdir(os.path.join(cfg.TEMPORAL, 'lumpy-sv'))
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
  assert common.execCmd('sudo cp bin/* /usr/local/bin/')[0]
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm -r lumpy-sv')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerLUMPY(cfg))

# Ensembl-VEP
def checkEVEP(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'ensembl-vep/vep'))
def checkVerEVEP(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'ensembl-vep/vep | grep ensembl-vep'), showcmd=False)
  assert ret[0], 'Ensembl-VEP is not installed.'
  ver = ret[1].split()[-1]
  return ver
def installEVEP(cfg, ver):
  print('Install Ensembl-VEP ...')
  assert common.execCmd("yes '' | sudo cpan Module::Build", verbose=True)[0], 'Module::Build install error.'
  assert common.execCmd("yes '' | sudo cpan Archive::Zip", verbose=True)[0], 'Archive::Zip install error.'
  assert common.execCmd("yes '' | sudo cpan DBD::mysql", verbose=True)[0], 'DBD::mysql install error.'
  assert common.execCmd("yes '' | sudo cpan DBI", verbose=True)[0], 'DBI install error.'
  assert common.execCmd("yes '' | sudo cpan Bio::Root::Version", verbose=True)[0], 'Bio::Root::Version install error.'
  os.chdir(cfg.APPS_DIR)  
  assert common.execCmd("git clone https://github.com/Ensembl/ensembl-vep.git")[0], 'Git clone error.'
  os.chdir(os.path.join(cfg.APPS_DIR, 'ensembl-vep'))
  assert common.execCmd("yes '' | perl INSTALL.pl", verbose=True)[0], 'Install error.'
  print('Completed.')
  print('>ver.', checkVerEVEP(cfg))

# SIFT-4G
def checkSift(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'SIFT4G_Annotator.jar'))
def checkVerSift(cfg):
  ret = common.execCmd('java -jar ' + os.path.join(cfg.APPS_DIR, 'SIFT4G_Annotator.jar -c 2>&1 >/dev/null | grep OPTIONS'), showcmd=False)
  ver = ret[1].split()[-2]
  return ver
def installSift(cfg, ver):
  print('Install SHIFT-4G...')
  os.chdir(cfg.APPS_DIR)
  assert common.execCmd('wget https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar')[0], 'Download error.'
  print('Completed.')
  print('>ver.', checkVerSift(cfg))
def installSiftDB(cfg, species, dbver):
  print('Install SHIFT-4G Database...')
  os.chdir(cfg.DB_DIR)
  assert common.execCmd('wget https://sift.bii.a-star.edu.sg/sift4g/public/' + species + '/' + dbver + '.zip')[0], 'Download error.'
  assert common.execCmd('unzip ' + dbver + '.zip'), 'Expansion error.'
  print('Completed.')
  
# HLA-HD
def checkHLAHD(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'hlahd/bin/hlahd.sh'))
def checkVerHLAHD(cfg):
  ret = common.execCmd('sh '+os.path.join(cfg.APPS_DIR, 'hlahd/bin/hlahd.sh 2>/dev/null | head -n 1'), showcmd=False)
  ver = ret[1].split()[-1].strip()
  return ver
def installHLAHD(cfg, ver):
  if not checkBowtie(cfg):
    print('Please install Bowtie2 before the installation of HLA-HD.')
  print('Install HLA-HD ...')
  if not os.path.exists(os.path.join(cfg.TEMPORAL, 'hlahd.'+ver+'.tar.gz')):
    print('HLA-HD source file was not found. Please download ' + 'hlahd.'+ver+'.tar.gz' + ' from "https://www.genome.med.kyoto-u.ac.jp/HLA-HD" and save it in $(TEMPORAL) directory.')
    return
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('tar -xvzf ' + os.path.join(cfg.TEMPORAL, 'hlahd.'+ver+'.tar.gz'))[0]
  os.chdir(os.path.join(cfg.TEMPORAL, 'hlahd.'+ver))
  assert common.execCmd('sh ./install.sh', verbose=True)[0], 'HLA-HD install error.'
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('mv ./hlahd.'+ver + ' ' + os.path.join(cfg.APPS_DIR, 'hlahd'), verbose=True)[0], 'HLA-HD install error.'
  print('Completed.')
  print('>ver.', checkVerHLAHD(cfg))

# InterProScan
def checkInterPro(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'interproscan/interproscan.sh'))
def checkVerInterPro(cfg):
  os.chdir(os.path.join(cfg.APPS_DIR, 'interproscan'))
  ret = common.execCmd('./interproscan.sh -version', showcmd=False)
  assert ret[0], 'InterProScan is not installed.'
  ver = ret[1].splitlines[0].replace('InterProScan version', '').strip()
  return ver
def installInterPro(cfg, ver):
  print('Install InterPro ...')
  assert common.execCmd('sudo apt-get install libdw1', verbose=True)[0], 'Install error.'
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/'+ver.split('.')[0]+'/'+ver+'/interproscan-'+ver+'-64-bit.tar.gz', verbose=True)[0], 'Download error.'
  assert common.execCmd('tar -xvzf interproscan-'+ver+'-64-bit.tar.gz', verbose=True)[0], 'Expansion error.'
  os.chdir(os.path.join(cfg.TEMPORAL, 'interproscan-'+ver))
  assert common.execCmd('python setup.py -f interproscan.properties', verbose=True)[0], 'Initialization error.'
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('mv interproscan-'+ver + os.path.join(cfg.APPS_DIR, 'interproscan'), verbose=True)[0], 'Initialization error.'
  assert common.execCmd('rm interproscan*', verbose=True)[0]
  print('Completed.')
  print('> ver.', checkVerInterPro(cfg))

# HTSeq
def checkHTSeq(cfg):
  ret = common.execCmd('which htseq-count', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerHTSeq(cfg):
  ret = common.execCmd('htseq-count --version', showcmd=False)
  assert ret[0], 'HTSeq is not installed.'
  ver = ret[1].strip()
  return ver
def installHTSeq(cfg, ver):
  print('Install HTSeq ...')
  assert common.execCmd('pip install HTSeq', verbose=True)[0], 'Install error.'
  print('Completed.')
  print('>ver.', checkVerHTSeq(cfg))
  
# featureCounts
def checkFeatCounts(cfg):
  ret = common.execCmd('which featureCounts', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerFeatCounts(cfg):
  ret = common.execCmd('featureCounts -v', showcmd=False)
  ver = ret[2].split()[1].strip()
  return ver
def installFeatCounts(cfg, ver):
  print('Install featureCounts...')
  assert common.execCmd('sudo apt-get update')[0]
  assert common.execCmd('sudo apt-get install -y subread', verbose=True)[0], 'Install error.'
  print('Completed.') 
  print('>ver.', checkVerFeatCounts(cfg))

# Cufflinks
def checkCuff(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR, 'cuff/cufflinks'))
def checkVerCuff(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR, 'cuff/cufflinks 2>&1 >/dev/null | head -n 1'), showcmd=False)
  ver = ret[1].split()[1]
  return ver
def installCuff(cfg, ver):
  print('Install Cufflinks...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-' + ver + '.Linux_x86_64.tar.gz')[0], 'Download error.'
  assert common.execCmd('tar xvf ./cufflinks-' + ver + '.Linux_x86_64.tar.gz')[0]
  os.makedirs(os.path.join(cfg.APPS_DIR, 'cuff'),exist_ok=True)
  assert common.execCmd('mv ./cufflinks-' + ver + '.Linux_x86_64/* '+cfg.APPS_DIR+'/cuff')[0]
  assert common.execCmd('rm -r ./cufflinks*')[0]
  os.chdir(cfg.WORK_SPACE)
  installCumme(cfg)
  print('Completed.')
  print('>ver.', checkVerCuff(cfg))

# RSEM
def checkRSEM(cfg):
  ret = common.execCmd('which rsem-calculate-expression', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerRSEM(cfg):
  ret = common.execCmd('rsem-calculate-expression --version', showcmd=False)
  assert ret[0], 'RSEM is not installed.'
  ver = ret[1].split()[-1]
  return ver
def installRSEM(cfg, ver):
  print('Install RSEM ...')
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('wget https://github.com/deweylab/RSEM/archive/v' + ver + '.tar.gz')[0], 'Download error.'
  assert common.execCmd('tar zxf ./v' + ver + '.tar.gz')[0], 'Expansion error.'
  os.chdir(os.path.join(cfg.TEMPORAL, 'RSEM-'+ver))
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
  assert common.execCmd('sudo make install', verbose=True)[0], 'Install error.'
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm v' + ver + '.tar.gz')[0]
  assert common.execCmd('rm -r RSEM*')[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerRSEM(cfg))

# BiocManager (R)
def checkBM(cfg):
  proc = subprocess.run('R --no-save --slave --vanilla < ' + os.path.join(cfg.SCRIPT_DIR, 'checkBiocManager.R'), stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0:
    return 'FALSE' not in proc.stdout
def installBM(cfg, ver):
  os.system('curl --output ' + os.path.join(cfg.SCRIPT_DIR, 'installBiocManager.R') + ' https://raw.githubusercontent.com/YujiSue/ysngs/main/R/installBiocManager.R')
  proc = subprocess.run('R --no-save --slave --vanilla < ' + os.path.join(cfg.SCRIPT_DIR, 'installBiocManager.R'), stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    print('> ',proc.stdout.splitlines()[-1])
# EdgeR (R)
def installEdgeR(cfg):
  os.system('curl --output ' + os.path.join(cfg.SCRIPT_DIR, 'installEdgeR.R') + ' https://raw.githubusercontent.com/YujiSue/ysngs/main/R/installEdgeR.R')
  proc = subprocess.run('R --no-save --slave --vanilla < ' + os.path.join(cfg.SCRIPT_DIR, 'installEdgeR.R'), stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    print('> ',proc.stdout.splitlines()[-1])
# CummeRbund (R)
def installCumme(cfg):
  os.system('curl --output ' + os.path.join(cfg.SCRIPT_DIR, 'installcummeR.R') + ' https://raw.githubusercontent.com/YujiSue/ysngs/main/R/installcummeR.R')
  proc = subprocess.run('R --no-save --slave --vanilla < ' + os.path.join(cfg.SCRIPT_DIR, 'installcummeR.R'), stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    print('> ',proc.stdout.splitlines()[-1])

# MACS2
def checkMACS(cfg):
  proc = subprocess.run('which macs2', stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    return os.path.exists(proc.stdout.splitlines()[0])
def checkVerMACS(cfg):
  ret = common.execCmd('macs2 --version | tail -n 1', showcmd=False)
  assert ret[0], 'MACS2 is not installed.'
  ver = ret[1].split()[-1]
  return ver
def installMACS(cfg, ver):
  print('Install MACS2 ...')
  os.chdir(cfg.TEMPORAL)
  print('  Downloads sources ...') 
  assert common.execCmd('wget https://github.com/macs3-project/MACS/archive/refs/tags/'+ver+'.tar.gz', verbose=True)[0], 'Download error.'
  assert common.execCmd('tar -xvzf '+ver+'.tar.gz', verbose=True)[0], 'Expansion error.'
  os.chdir(os.path.join(cfg.TEMPORAL,'MACS-'+ver[1:]))
  assert common.execCmd("sed -i 's/install_requires = \[f\"numpy>={numpy_requires}\",\]/install_requires = \[f\"numpy{numpy_requires}\",\]/' ./setup.py", verbose=True)[0]
  assert common.execCmd("pip install -e .", verbose=True)[0]
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm '+ver+'.tar.gz', verbose=True)[0]
  assert common.execCmd('rm -r MACS*', verbose=True)[0]  
  print('Completed.')
  print('>ver.', checkVerMACS(cfg))

# MEME
def checkMEME(cfg):
  return os.path.exists(os.path.join(cfg.APPS_DIR,'meme/bin/meme'))
def checkVerMEME(cfg):
  ret = common.execCmd(os.path.join(cfg.APPS_DIR,'meme/bin/meme -version'), showcmd=False)
  assert ret[0], 'MEME Suite is not installed.'
  ver = ret[1].strip()
  return ver
def installMEME(cfg, ver):
  print('Install MEME Suite ...')
  os.chdir(cfg.TEMPORAL)
  print('  Downloads sources ...') 
  assert common.execCmd('wget https://meme-suite.org/meme/meme-software/'+ver+'/meme-'+ver+'.tar.gz', verbose=True)[0]
  assert common.execCmd('tar -xvzf ./meme-'+ver+'.tar.gz', verbose=True)[0]
  os.chdir('meme-'+ver)
  assert common.execCmd('sudo cpan force install Sys::Info',verbose=True)[0]
  assert common.execCmd('./configure --prefix='+os.path.join(cfg.APPS_DIR, 'meme --enable-build-libxml2 --enable-build-libxslt'), verbose=True)[0]
  assert common.execCmd('make', verbose=True)[0]
  #assert common.execCmd('make test', verbose=True)[0]
  assert common.execCmd('sudo make install', verbose=True)[0]
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('rm -r ./meme*', verbose=True)[0]
  os.chdir(cfg.WORK_SPACE)
  print('Completed.')
  print('>ver.', checkVerMEME(cfg))

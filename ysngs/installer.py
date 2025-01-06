import os
import subprocess
from subprocess import PIPE 
from ysngs import common

# Moirei
def checkMoirei():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'bin', 'moirei'))
def checkVerMoirei():
  res = common.execCmd(f"{os.path.join(os.environ['HYM_APP'], 'bin', 'moirei')} --version", showcmd=False)
  return res[1]
def installMoirei(prop):
  if checkMoirei():
    print('moirei is installed.')
  else:
    print('Install moirei ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.gitClone(prop['url'])
    os.chdir('Moirei')
    assert common.execCmd(f"cmake -DINSTALL_SLIB=ON -S . -B build", showcmd=False, verbose=True)[0], "Failed to cmake configure."
    assert common.execCmd(f"cmake --build build", showcmd=False, verbose=True)[0], "Failed to build."
    assert common.execCmd(f"cmake --install build --prefix {os.environ['HYM_APP']}", showcmd=False, verbose=True)[0], "Failed to build."
    print('Completed.')
    print('> ver.', checkVerMoirei())

# Sutoku
def checkSutoku():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'bin', 'sutoku'))
def checkVerSutoku():
  res = common.execCmd(f"{os.path.join(os.environ['HYM_APP'], 'bin', 'sutoku')} --version", showcmd=False)
  return res[1]
def installSutoku(prop):
  if checkSutoku():
    print('sutoku is installed.')
  else:
    print('Install sutoku ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.gitClone(prop['url'])
    os.chdir('Sutoku')
    assert common.execCmd(f"cmake -DINSTALL_SLIB=ON -S . -B build", showcmd=False, verbose=True)[0], "Failed to cmake configure."
    assert common.execCmd(f"cmake --build build", showcmd=False, verbose=True)[0], "Failed to build."
    assert common.execCmd(f"cmake --install build --prefix {os.environ['HYM_APP']}", showcmd=False, verbose=True)[0], "Failed to build."
    print('Completed.')
    print('> ver.', checkVerSutoku())
    
# Cromwell
def checkCromwell():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'cromwell.jar'))
def checkWom():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'womtool.jar'))
def checkVerCromwell():
  res = common.execCmd(f"java -jar {os.path.join(os.environ['HYM_APP'], 'cromwell.jar')} --version", showcmd=False)
  #return res[1].split(' ')[-1]
  return res[1]
def checkVerWom():
  res = common.execCmd(f"java -jar {os.path.join(os.environ['HYM_APP'], 'womtool.jar')} --version", showcmd=False)
  #return res[1].split(' ')[-1]
  return res[1]
def installCromwell(prop):
  if checkCromwell():
    print('Cromwell is installed.')
  else:
    print('Install cromwell ...')
    common.curlDownload(prop['url'][0], output=os.path.join(os.environ['HYM_APP'], 'cromwell.jar'))
    print('Completed.')
    print('> ver.', checkVerCromwell())
  if checkWom():
    print('Womtool is installed.')
  else:
    print('Install womtool ...')
    common.curlDownload(prop['url'][1], output=os.path.join(os.environ['HYM_APP'], 'womtool.jar'))
    print('Completed.')
    print('> ver.', checkVerWom())

# Miniconda
def checkConda():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'miniconda', 'bin', 'conda'))
def checkVerConda():
  res = common.execCmd(f"{os.path.join(os.environ['HYM_APP'], 'miniconda', 'bin', 'conda')} --version", showcmd=False)
  return res[1].split(' ')[-1]
def installConda(prop):
  if checkConda():
    print('Conda is already installed.')
  else:
    print('Install miniconda ...')
    os.chdir(os.environ['HYM_APP'])
    os.makedirs('miniconda', exist_ok=True)
    common.curlDownload(prop['url'], output=os.path.join('miniconda', 'install.sh'))
    assert common.execCmd(f"bash miniconda/install.sh -b -u -p miniconda", showcmd=False, verbose=True)[0], "Failed to run install shell script."
    print('Completed.')
    print('> ver.', checkVerConda())

# IGV
def checkIGV():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'igv/igv.sh')) 
def checkVerIGV():
  res = common.execCmd(os.path.join(os.environ['HYM_APP'], 'igv/igv.sh --version'), showcmd=False)
  assert res[0], 'IGV is not installed.'
  return res[1].splitlines()[-1]
def installIGV(prop):
  if checkIGV():
    print('IGV is already installed.')
  else:
    print('Install IGV...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[1]
    assert common.execCmd(f"unzip {fname}", showcmd=False)[0], 'Expansion error.'
    os.makedirs(os.path.join(os.environ['HYM_APP'], 'igv'), exist_ok=True)
    assert common.execCmd(f"mv IGV_Linux* {os.path.join(os.environ['HYM_APP'], 'igv')}")[0], 'File move error.'
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('>ver.', checkVerIGV())

# SRA-Toolkit
def checkSRA():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'sra', 'test-sra'))
def checkVerSRA():
  res = common.execCmd(f"{os.path.join(os.environ['HYM_APP'], 'sra', 'test-sra')} | grep 'NCBI SRA Toolkit release version:' | head -n 1", showcmd=False)
  assert res[0], 'SRA-Toolkit is not installed.'
  return res[1].replace('NCBI SRA Toolkit release version:','').replace('.<br/>','')
def installSRA(prop):
  if checkSRA():
    print('SRA-Toolkit installed.')
  else:
    print('Install SRA-Toolkit...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[1]
    assert common.execCmd(f"tar xvf {fname}", showcmd=False)[0]
    os.makedirs(os.path.join(os.environ['HYM_APP'], 'sra'), exist_ok=True)
    assert common.execCmd(f"cp -r ./sratoolkit*/bin/* {os.path.join(os.environ['HYM_APP'], 'sra')}", showcmd=False)[0]
    assert common.execCmd('rm -r ./sratoolkit*', showcmd=False)[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('> ver.', checkVerSRA())
    common.execCmd(f"{os.path.join(os.environ['HYM_APP'], 'sra', 'vdb-config')} --interactive")

# Samtools
def checkSamtools():
  res = common.execCmd('which samtools', showcmd=False)
  return res[0] and os.path.exists(res[1])
def checkVerSamtools():
  res = common.execCmd('samtools --version | head -n 1', showcmd=False)
  assert res[0], 'Samtools is not installed.'
  return res[1].replace('samtools','')
def installSamtools(prop):
  if checkSamtools():
    print('samtools is installed.')
  else:
    print('Install Samtools...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[1]
    assert common.execCmd(f"tar xvf {fname}", showcmd = False)[0], 'Expansion error.'
    os.chdir(fname[0:fname.find('.tar')])
    assert common.execCmd('./configure', showcmd = False, verbose=True)[0], 'Configure error.'
    assert common.execCmd('make -j8', showcmd = False, verbose=True)[0], 'Make error.'
    assert common.execCmd('sudo make install', showcmd = False, verbose=True)[0], 'Install error.'
    print('Completed.')
    print('> ver.', checkVerSamtools())

# BCFtools
def checkBCFtools():
  res = common.execCmd('which bcftools', showcmd=False)
  return res[0] and os.path.exists(res[1])
def checkVerBCFtools():
  res = common.execCmd('bcftools --version | head -n 1', showcmd=False)
  assert res[0], 'BCFtools is not installed.'
  return res[1].replace('bcftools','')
def installBCFtools(prop):
  if checkBCFtools():
    print('BCFtools is installed.')
  else:
    print('Install BCFtools ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[1]
    assert common.execCmd(f"tar xvf {fname}", showcmd=False)[0], 'Expansion error.'
    os.chdir(fname[0:fname.find('.tar')])
    assert common.execCmd('./configure', showcmd=False, verbose=True)[0], 'Configure error.'
    assert common.execCmd('make -j8', showcmd=False, verbose=True)[0], 'Make error.'
    assert common.execCmd('sudo make install', showcmd=False)[0], 'Install error.'
    os.chdir(os.environ['HYM_TEMP'])
    assert common.execCmd('rm -r ./bcftools*', showcmd=False)[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('> ver.', checkVerBCFtools())

# Picard
def checkPicard():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'picard.jar'))
def checkVerPicard():
  res = common.execCmd('java -jar ' + os.path.join(os.environ['HYM_APP'], 'picard.jar MarkDuplicates --version'), showcmd=False)
  return res[2].replace('Version:','')
def installPicard(prop):
  if checkPicard():
    print('Picard is installed.')
  else:
    print('Install Picard ...')  
    common.curlDownload(prop['url'], output=os.path.join(os.environ['HYM_APP'], 'picard.jar'))
    print('Completed.')
    os.chdir(os.environ['HYM_WS'])
    print('>ver.', checkVerPicard())
  
# GATK
def checkGATK():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'gatk'))
def checkVerGATK():
  res = common.execCmd(os.path.join(os.environ['HYM_APP'], 'gatk/gatk --list'), showcmd=False)
  return res[2].splitlines()[0].split('-')[-2]
def installGATK(prop):
  if checkSamtools():
    print('GATK is installed.')
  else:
    print('Install GATK ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[1]
    assert common.execCmd(f"unzip -o {fname}", showcmd=False)[0], 'Expansion error.'
    assert common.execCmd(f"mv ./gatk* {os.path.join(os.environ['HYM_APP'], 'gatk')}", showcmd=False)[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('>ver.', checkVerGATK())
  
# CutAdapt
def checkCut():
  res = common.execCmd('which cutadapt', showcmd=False)
  return res[0] and os.path.exists(res[1].strip())
def checkVerCut():
  res = common.execCmd('cutadapt --version', showcmd=False)
  ver = res[1].strip()
  return ver
def installCut(prop):
  if checkCut():
    print('cutadapt is installed.')
  else:
    print('Install cutadapt  ...')
    assert common.execCmd('pip install cutadapt', verbose=True)[0], 'Install error.'
    print('Completed.')
    print('>ver.', checkVerCut())

# FastQC
def checkFQC():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'FastQC/fastqc'))
def checkVerFQC():
  res = common.execCmd(os.path.join(os.environ['HYM_APP'], 'FastQC/fastqc --version &'), showcmd=False)
  assert res[0], 'fastqc is not installed.'
  return res[1].replace('FastQC ','').strip()
def installFQC(prop):
  if checkFQC():
    print('FastQC is installed.')
  else:
    print('Install FastQC  ...')
    os.chdir(os.environ['HYM_TEMP'])
    assert common.curlDownload(prop['url'])[0], 'Download error.'
    fname = os.path.split(prop['url'])[1]
    assert common.execCmd(f"unzip {fname}")[0], 'Expansion error.'
    assert common.execCmd('mv FastQC $HYM_APP')[0], 'File move error.'
    assert common.execCmd('chmod a+x $HYM_APP/FastQC/fastqc')[0], 'Permission setting error.'
    assert common.execCmd('rm -r ./*')[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('>ver.', checkVerFQC())

# fastp
def checkFP():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'fastp'))
def checkVerFP():
  res = common.execCmd(os.path.join(os.environ['HYM_APP'], 'fastp --version'), showcmd=False) 
  ver = res[2].replace('fastp','').strip()
  return ver
def installFP(prop):
  if checkFP():
    print('fastp is installed.')
  else:
    print('Install fastp  ...')
    assert common.curlDownload(prop['url'], '$HYM_APP/fastp')[0], 'Download error.'
    assert common.execCmd('chmod a+x $HYM_APP/fastp')[0], 'Permission setting error.'
    print('Completed.')
    print('>ver.', checkVerFP())

# BWA
def checkBWA():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'bwa'))
def checkVerBWA():
  res = common.execCmd(os.path.join(os.environ['HYM_APP'], 'bwa 2>&1 >/dev/null | grep Version:'), showcmd=False)
  return res[1].replace('Version:', '')
def installBWA(prop):
  if checkBWA():
    print('bwa is installed.')
  else:
    print('Install BWA  ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[-1]
    assert common.execCmd(f"tar xvf {fname}", showcmd=False)[0], 'Expansion error.'
    fname = f"bwa-{fname[1:fname.find('.tar')]}"
    os.chdir(fname)
    common.execCmd(f"sed -i 's/const uint8_t rle_auxtab/extern const uint8_t rle_auxtab/g' rle.h", showcmd=False)
    assert common.execCmd('make -j8', showcmd=False, verbose=True)[0], 'Make error.'
    assert common.execCmd('cp bwa $HYM_APP', showcmd=False)[0], 'File copy error.'
    assert common.execCmd(f"rm -r {os.path.join(os.environ['HYM_TEMP'], '*')}", showcmd=False)[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('> ver.', checkVerBWA())

# Bowtie2
def checkBowtie():
  return os.path.exists(os.path.join(os.environ['HYM_APP'],'bowtie2','bowtie2'))
def checkVerBowtie():
  res = common.execCmd('bowtie2 --version | head -n 1', showcmd=False)
  assert res[0], 'Bowtie2 is not installed.'
  return res[1].split()[-1]
def installBowtie(prop):
  if checkBowtie():
    print('Bowtie2 is installed.')
  else:
    print('Install Bowtie2 ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[-1]
    assert common.execCmd(f"unzip -o {fname}", showcmd=False)[0], 'Expansion error.'
    fname = fname[0:fname.find('.zip')]
    assert common.execCmd(f"mv {fname} $HYM_APP/bowtie2", showcmd=False, verbose=True)[0], 'Make error.'
    assert common.execCmd('rm -r ./*', showcmd=False)[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('>ver.', checkVerBowtie())

# STAR
def checkSTAR():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'STAR'))
def checkVerSTAR():
  res = common.execCmd(os.path.join(os.environ['HYM_APP'], 'STAR | grep version'), showcmd=False)
  return res[1].split('=')[-1]
def installSTAR(prop):
  if checkSTAR():
    print('STAR is installed.')
  else:
    print('Install STAR ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[-1]
    assert common.execCmd(f"tar xvf {fname}", showcmd=False)[0], 'Expansion error.'
    fname = f"STAR-{fname[0:fname.find('.tar')]}"
    os.chdir(os.path.join(fname, 'source'))
    assert common.execCmd('make -j8 STAR', showcmd=False, verbose=True)[0], 'Make error.'
    os.chdir(os.environ['HYM_TEMP'])
    assert common.execCmd(f"mv {fname} {os.path.join(os.environ['HYM_APP'], 'STAR')}")[0]
    assert common.execCmd('rm -r ./*', showcmd=False)[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('>ver.', checkVerSTAR())
  
# HISAT2
def checkHISAT2():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'hisat2', 'hisat2'))
def checkVerHISAT2():
  res = common.execCmd(os.path.join(os.environ['HYM_APP'], 'hisat2', 'hisat2 --version | head -n 1'), showcmd=False)
  return res[1].split(' ')[-1]
def installHISAT2(prop):
  if checkHISAT2():
    print('HISAT2 is installed.')
  else:
    print('Install HISAT2...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'], 'hisat2.zip')    
    assert common.execCmd('unzip hisat2.zip')[0], 'Expansion error.'
    assert common.execCmd('rm hisat2.zip')[0]
    assert common.execCmd(f"mv hisat2* {os.path.join(os.environ['HYM_APP'], 'hisat2')}")[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('>ver.', checkVerHISAT2())

# IVC
def checkStrelka():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'strelka'))
def checkVerStrelka():
  ret = common.execCmd('python2 '+os.path.join(os.environ['HYM_APP'], 'strelka/bin/configureStrelkaGermlineWorkflow.py --version'), showcmd=False)
  assert ret[0], 'Strelka is not installed.'
  ver = ret[1].strip()
  return ver
def installStrelka(prop):
  print('Install Strelka...')
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('wget https://github.com/Illumina/strelka/releases/download/v' + ver + '/strelka-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Download error.'
  assert common.execCmd('tar xvjf strelka-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Expansion error.'
  assert common.execCmd('mv -f strelka-' + ver + '.centos6_x86_64' + ' ' + os.path.join(os.environ['HYM_APP'], 'strelka'))[0]
  assert common.execCmd('rm -r strelka*')[0]
  os.chdir(os.environ['HYM_WS'])
  print('Completed.')
  print('>ver.', checkVerStrelka())
 
# Manta
def checkManta(cfg):
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'Manta/bin/configManta.py'))
def checkVerManta(cfg):
  ret = common.execCmd(os.path.join(os.environ['HYM_APP'], 'Manta/bin/configManta.py --version'), showcmd=False)
  assert ret[0], 'Manta is not installed.'
  return ret[1].strip()
def installManta(cfg, ver):
  print('Install Manta...')
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('wget https://github.com/Illumina/manta/releases/download/v' + ver + '/manta-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Download error.'
  assert common.execCmd('tar xvjf manta-' + ver + '.centos6_x86_64.tar.bz2')[0], 'Expansion error.'
  assert common.execCmd('mv ./manta-' + ver + '.centos6_x86_64 ' + os.path.join(os.environ['HYM_APP'], 'Manta'))[0]
  assert common.execCmd('rm -r ./manta*')[0]
  os.chdir(os.environ['HYM_WS'])
  print('Completed.')
  print('>ver.', checkVerManta(cfg))

# TVC
def checkTVC():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'TVC/bin/tvc'))
def checkVerTVC():
  ret = common.execCmd(os.path.join(os.environ['HYM_APP'], 'TVC/bin/tvc') + ' --version')
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
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('wget http://updates.iontorrent.com/tvc_standalone/tvc-' + ver+'.tar.gz')[0], 'Download error.'
  build = os.path.join(os.environ['HYM_TEMP'], 'temp')
  os.makedirs(build, exist_ok=True)
  assert common.execCmd('cp '+os.environ['HYM_TEMP']+'/tvc-' + ver + '.tar.gz ' + build + '/tvc-' + ver + '.tar.gz')[0]
  install = os.path.join(os.environ['HYM_APP'], 'TVC')
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
  assert common.execCmd(os.environ['HYM_APP']+'/TVC/bin/samtools bam2fq '+cfg.TEST_DIR+'/test.bam > '+cfg.TEST_DIR+'/test.fq')[0]
  assert common.execCmd('rm -r '+build)[0]
  assert common.execCmd('rm '+os.environ['HYM_TEMP']+'/tvc*')[0]
  os.chdir(os.environ['HYM_WS'])
  print('Completed.')
  print('>ver.', checkVerTVC(cfg))

# GDV
def checkGDV():
  res = common.execCmd('docker images')
  return res[0] and 'google/deepvariant' in res[1]
def checkVerGDV():
  res = common.execCmd('docker images | grep google/deepvariant')
  assert res[0], 'Google DeepVariant is not installed.'
  return res[1].split()[1]
def installGDV(prop):
  if checkGDV():
    print('DeepVariant is installed.')
  else:
    print('Install DeepVariant  ...')
    assert common.execCmd(f"sudo docker pull google/deepvariant:{prop['ver']}")[0]
    print('Completed.')
    print('>ver.', checkVerGDV())

# CNVnator
def checkCNVnator(cfg):
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'CNVnator', 'cnvnator'))
def checkVerCNVnator(cfg):
  ret = common.execCmd(os.path.join(os.environ['HYM_APP'], 'CNVnator', 'cnvnator 2>&1 >/dev/null | grep CNVnator'), showcmd=False)
  assert ret[0], 'CNVnator is not installed.'
  ver = ret[1].split()[1].strip()
  return ver
def installCNVnator(cfg, ver):
  print('Install CNVnator...')
  assert common.execCmd('echo deb http://security.ubuntu.com/ubuntu xenial-security main | sudo tee -a /etc/apt/sources.list')[0]
  assert common.execCmd('sudo apt-get update')[0]
  assert common.execCmd('sudo apt-get install libssl1.0.0 libreadline-dev')[0]
  os.chdir(os.environ['HYM_APP'])
  assert common.execCmd('git clone --recursive https://github.com/abyzovlab/CNVnator.git')[0], 'Download error.'
  os.chdir(os.path.join(os.environ['HYM_APP'], 'CNVnator'))
  assert common.download('https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2')[0]
  assert common.execCmd('tar -xvf samtools-1.2.tar.bz2')[0]
  assert common.execCmd('rm samtools-1.2.tar.bz2')[0]
  os.chdir(os.path.join(os.environ['HYM_APP'], 'CNVnator', 'samtools-1.2'))
  assert common.execCmd('make -j8')[0], 'Samtools make failed.'
  os.chdir(os.path.join(os.environ['HYM_APP'], 'CNVnator'))
  assert common.execCmd('ln -s samtools-1.2 samtools')[0]
  os.environ['ROOTSYS'] = os.path.join(os.environ['HYM_APP'], 'root')
  assert common.execCmd('make -j8')[0], 'Make failed.'
  os.chdir(os.path.join(os.environ['HYM_WS']))
  print('Completed.')
  print('>ver.', checkVerCNVnator(cfg))

# DELLY
def checkDELLY(cfg):
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'delly/src/delly'))
def checkVerDELLY(cfg):
  ret = common.execCmd(os.path.join(os.environ['HYM_APP'] ,'delly/src/delly 2>&1 >/dev/null | grep Version:'), showcmd=False)
  ver = ret[1].split()[-1][:-1]
  return ver
def installDELLY(cfg, ver):
  print('Install DELLY...')
  os.chdir(os.environ['HYM_APP'])
  assert common.execCmd('git clone --recursive https://github.com/dellytools/delly.git')[0], 'Download error.'
  os.chdir('./delly')
  assert common.execCmd('make all -j8')[0], 'Compile error.'
  os.chdir(os.environ['HYM_WS'])
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
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('git clone --recursive https://github.com/arq5x/lumpy-sv.git')[0], 'Download error.'
  os.chdir(os.path.join(os.environ['HYM_TEMP'], 'lumpy-sv'))
  assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
  assert common.execCmd('sudo cp bin/* /usr/local/bin/')[0]
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('rm -r lumpy-sv')[0]
  os.chdir(os.environ['HYM_WS'])
  print('Completed.')
  print('>ver.', checkVerLUMPY(cfg))

# Ensembl-VEP
def checkEVEP():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'ensembl-vep/vep'))
def checkVerEVEP():
  ret = common.execCmd(os.path.join(os.environ['HYM_APP'], 'ensembl-vep/vep | grep ensembl-vep'), showcmd=False)
  assert ret[0], 'Ensembl-VEP is not installed.'
  return ret[1].split()[-1]
def installEVEP(prop):
  if checkEVEP():
    print('Ensembl-VEP is installed.')
  else:
    print('Install Ensembl-VEP ...')
    assert common.execCmd("yes '' | sudo cpan Module::Build", showcmd=False, verbose=True)[0], 'Module::Build install error.'
    assert common.execCmd("yes '' | sudo cpan Archive::Zip", showcmd=False, verbose=True)[0], 'Archive::Zip install error.'
    assert common.execCmd("yes '' | sudo cpan DBD::mysql", showcmd=False, verbose=True)[0], 'DBD::mysql install error.'
    assert common.execCmd("yes '' | sudo cpan DBI", showcmd=False, verbose=True)[0], 'DBI install error.'
    assert common.execCmd("yes '' | sudo cpan Bio::Root::Version", showcmd=False, verbose=True)[0], 'Bio::Root::Version install error.'
    os.chdir(os.environ['HYM_APP'])  
    assert common.execCmd(f"git clone {prop['url']}", showcmd=False)[0], 'Git clone error.'
    os.chdir(os.path.join(os.environ['HYM_APP'], 'ensembl-vep'))
    assert common.execCmd("yes '' | perl INSTALL.pl", showcmd=False, verbose=True)[0], 'Install error.'
    print('Completed.')
    print('>ver.', checkVerEVEP())

# SIFT-4G
def checkSift():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'SIFT4G_Annotator.jar'))
def checkVerSift():
  ret = common.execCmd('java -jar $HYM_APP/SIFT4G_Annotator.jar -c 2>&1 >/dev/null | grep OPTIONS', showcmd=False)
  ver = ret[1].split()[-2]
  return ver
def installSift(prop):
  if checkSift():
    print('SIFT-4G is installed.')
  else:
    print('Install SIFT-4G ...')
    common.curlDownload(prop['url'], '$HYM_APP')
    print('Completed.')
    print('>ver.', checkVerSift())
  
# HLA-HD
def checkHLAHD():
  return os.path.exists(os.path.join(os.environ['HYM_APP'],'hlahd/bin/hlahd.sh'))
def checkVerHLAHD():
  ret = common.execCmd('sh $HYM_APP/hlahd/bin/hlahd.sh 2>/dev/null | head -n 1', showcmd=False)
  ver = ret[1].split()[-1].strip()
  return ver
def installHLAHD(prop):
  if not checkBowtie():
    print('Please install Bowtie2 before the installation of HLA-HD.')
  print('Install HLA-HD ...')
  if not os.path.exists(os.path.join(os.environ['HYM_TEMP'], 'hlahd.'+ver+'.tar.gz')):
    print('HLA-HD source file was not found. Please download ' + 'hlahd.'+ver+'.tar.gz' + ' from "https://www.genome.med.kyoto-u.ac.jp/HLA-HD" and save it in $(TEMPORAL) directory.')
    return
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('tar -xvzf ' + os.path.join(os.environ['HYM_TEMP'], 'hlahd.'+ver+'.tar.gz'))[0]
  os.chdir(os.path.join(os.environ['HYM_TEMP'], 'hlahd.'+ver))
  assert common.execCmd('sh ./install.sh', verbose=True)[0], 'HLA-HD install error.'
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('mv ./hlahd.'+ver + ' ' + os.path.join(os.environ['HYM_APP'], 'hlahd'), verbose=True)[0], 'HLA-HD install error.'
  print('Completed.')
  print('>ver.', checkVerHLAHD(cfg))

# InterProScan
def checkInterPro():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'interproscan/interproscan.sh'))
def checkVerInterPro():
  os.chdir(os.path.join(os.environ['HYM_APP'], 'interproscan'))
  res = common.execCmd('./interproscan.sh -version', showcmd=False)
  assert res[0], 'InterProScan is not installed.'
  return res[1].splitlines[0].replace('InterProScan version', '')
def installInterPro(prop):
  if checkInterPro():
    print('InterProScan is installed.')
  else:
    print('Install InterProScan ...')
    assert common.execCmd('sudo apt-get install libdw1', showcmd=False, verbose=True)[0], 'Install error.'
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[-1]
    assert common.execCmd(f"tar xvf {fname}", showcmd=False)[0], 'Expansion error.'
    fname = fname[0:fname.find('.tar')]
    os.chdir(fname)
    assert common.execCmd('python setup.py -f interproscan.properties', showcmd=False, verbose=True)[0], 'Initialization error.'
    os.chdir(os.environ['HYM_TEMP'])
    assert common.execCmd(f"mv {fname} {os.path.join(os.environ['HYM_APP'], 'interproscan')}", showcmd=False, verbose=True)[0], 'Initialization error.'
    assert common.execCmd(f"rm interproscan*", showcmd=False, verbose=True)[0]
    print('Completed.')
    print('> ver.', checkVerInterPro())
  
# HTSeq
def checkHTSeq():
  ret = common.execCmd('which htseq-count', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerHTSeq():
  ret = common.execCmd('htseq-count --version', showcmd=False)
  assert ret[0], 'HTSeq is not installed.'
  ver = ret[1].strip()
  return ver
def installHTSeq(prop):
  print('Install HTSeq ...')
  assert common.execCmd('pip install HTSeq', verbose=True)[0], 'Install error.'
  print('Completed.')
  print('>ver.', checkVerHTSeq())
  
# featureCounts
def checkFeatCounts():
  res = common.execCmd('which featureCounts', showcmd=False)
  return res[0] and os.path.exists(res[1].strip())
def checkVerFeatCounts():
  ret = common.execCmd('featureCounts -v', showcmd=False)
  ver = ret[2].split()[1].strip()
  return ver
def installFeatCounts(prop):
  if checkFeatCounts():
    print('featureCounts is installed.')
  else:
    print('Install featureCounts...')
    assert common.execCmd('sudo apt-get update', showcmd=False)[0]
    assert common.execCmd('sudo apt-get install -y subread', showcmd=False, verbose=True)[0], 'Install error.'
    print('Completed.') 
    print('>ver.', checkVerFeatCounts())

# Cufflinks
def checkCuff():
  return os.path.exists(os.path.join(os.environ['HYM_APP'], 'cuff/cufflinks'))
def checkVerCuff():
  ret = common.execCmd(os.path.join(os.environ['HYM_APP'], 'cuff/cufflinks 2>&1 >/dev/null | head -n 1'), showcmd=False)
  ver = ret[1].split()[1]
  return ver
def installCuff(prop):
  if checkCuff():
    print('Cufflinks is installed.')
  else:
    print('Install Cufflinks...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[-1]
    assert common.execCmd(f"tar xvf {fname}", showcmd=False)[0], 'Expansion error.'
    fname = fname[0:fname.find('.tar')]
    assert common.execCmd(f"mv {fname} {os.path.join(os.environ['HYM_APP'],'cuff')}")[0]
    assert common.execCmd('rm -r ./cuff*')[0]
    os.chdir(os.environ['HYM_WS'])
    installCumme(prop)
    print('Completed.')
    print('>ver.', checkVerCuff())

# RSEM
def checkRSEM():
  ret = common.execCmd('which rsem-calculate-expression', showcmd=False)
  return ret[0] and os.path.exists(ret[1].strip())
def checkVerRSEM():
  ret = common.execCmd('rsem-calculate-expression --version', showcmd=False)
  assert ret[0], 'RSEM is not installed.'
  ver = ret[1].split()[-1]
  return ver
def installRSEM(prop):
  if checkRSEM():
    print('RSEM is installed.')
  else:  
    print('Install RSEM ...')
    os.chdir(os.environ['HYM_TEMP'])
    common.curlDownload(prop['url'])
    fname = os.path.split(prop['url'])[-1]
    assert common.execCmd(f"tar xvf {fname}", showcmd=False)[0], 'Expansion error.'
    fname = fname[0:fname.find('.tar')].replace('v', 'RSEM-')
    os.chdir(fname)
    assert common.execCmd('make -j8', verbose=True)[0], 'Make error.'
    assert common.execCmd('sudo make install', verbose=True)[0], 'Install error.'
    os.chdir(os.environ['HYM_TEMP'])
    assert common.execCmd('rm -r ./*')[0]
    os.chdir(os.environ['HYM_WS'])
    print('Completed.')
    print('>ver.', checkVerRSEM())

# BiocManager (R)
def checkBM():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkBiocManager.R')
#  if not os.path.exists(scrpt):
#    common.download('https://raw.githubusercontent.com/YujiSue/ysngs/main/R/checkBiocManager.R', output = os.path.join(cfg.SCRIPT_DIR, 'checkBiocManager.R'))
  ret = common.runRScript(scrpt, output = None, args = [])
  return ret[0] and 'FALSE' not in ret[1]
def checkVerBM():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkBiocManagerVer.R')
  #if not os.path.exists(os.path.join(cfg.SCRIPT_DIR, 'checkBiocManagerVer.R')):
  #  common.download('https://raw.githubusercontent.com/YujiSue/ysngs/main/R/checkBiocManagerVer.R', output = os.path.join(cfg.SCRIPT_DIR, 'checkBiocManagerVer.R'))
  ret = common.runRScript(scrpt, output = None, args = [])
  return ret[1].strip().split(' ')[-1][1:-1]
def installBM(prop):
  print('Install BiocManager(R) ...')
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installBiocManager.R')
  
  #if not os.path.exists(os.path.join(cfg.SCRIPT_DIR, 'installBiocManager.R')):
  #  common.download('https://raw.githubusercontent.com/YujiSue/ysngs/main/R/installBiocManager.R', output = os.path.join(cfg.SCRIPT_DIR, 'installBiocManager.R'))
  assert common.runRScript(scrpt, output = None, args = [])[0]
  print('Completed.')
  print('>ver.', checkVerBM())

# DESeq2 (R)
def checkDESeq():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkBMPkg.R')
  ret = common.runRScript(os.path.join(os.environ['HYM_SCRIPT'], 'checkBMPkg.R'), args=['DESeq2'], output = None)
  return ret[0] and 'TRUE' in ret[1]
def checkVerDESeq():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkPkgVer.R')
  ret = common.runRScript(os.path.join(os.environ['HYM_SCRIPT'], 'checkBMPkgVer.R'), args=['DESeq2'], output = None)
  return ret[1].strip()[1:-1]
def installDESeq(prop):
  print('Install DESeq2 (R) ...')
  assert common.runRScript(os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installBMPkg.R'), args=['DESeq2'], output = None)
  print('Completed.')
  print('>ver.', checkVerDESeq())

# EdgeR (R)
def checkEdgeR():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkBMPkg.R')
  ret = common.runRScript(os.path.join(os.environ['HYM_SCRIPT'], 'checkBMPkg.R'), args=['edgeR'], output = None)
  return ret[0] and 'TRUE' in ret[1]
def checkVerEdgeR():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkPkgVer.R')
  ret = common.runRScript(os.path.join(os.environ['HYM_SCRIPT'], 'checkPkgVer.R'), args=['edgeR'], output = None)
  return ret[1].strip().split(' ')[-1][1:-1]
def installEdgeR(prop):
  print('Install EdgeR (R) ...')
  assert common.runRScript(os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installBMPkg.R'), args=['edgeR'], output = None)
  print('Completed.')
  print('>ver.', checkVerEdgeR())

# CummeRbund (R)
def checkCumme():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkBMPkg.R')
#  if not os.path.exists(os.path.join(cfg.SCRIPT_DIR, 'checkBMPkg.R')):
#    common.download('https://raw.githubusercontent.com/YujiSue/ysngs/main/R/checkBMPkg.R', output = os.path.join(cfg.SCRIPT_DIR, 'checkBMPkg.R'))
  ret = common.runRScript(scrpt, args=['cummeRbund'], output = None)
  return ret[0] and 'TRUE' in ret[1]
def checkVerCumme():
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'checkPkgVer.R')
#  if not os.path.exists(os.path.join(cfg.SCRIPT_DIR, 'checkPkgVer.R')):
#    common.download('https://raw.githubusercontent.com/YujiSue/ysngs/main/R/checkPkgVer.R', output = os.path.join(cfg.SCRIPT_DIR, 'checkPkgVer.R'))
  ret = common.runRScript(scrpt, args=['cummeRbund'], output = None)
  return ret[1].strip().split(' ')[-1][1:-1]
def installCumme(prop):
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installBMPkg.R')
  print('Install CummeRbund (R) ...')
#  if not os.path.exists(os.path.join(cfg.SCRIPT_DIR, 'installBMPkg.R')):
#    common.download('https://raw.githubusercontent.com/YujiSue/ysngs/main/R/installBMPkg.R', output = os.path.join(cfg.SCRIPT_DIR, 'installBMPkg.R'))
  assert common.runRScript(scrpt, args=['cummeRbund', 'ggplot2', 'RSQLite'], output = None)
  print('Completed.')
  print('>ver.', checkVerCumme())

# Org DB
def installRHumanGeneSet(prop):
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installBMPkg.R')
  print('Install clusterProfiler (R) ...')
  assert common.runRScript(scrpt, args=['org.Hs.eg.db'], output = None)[0]
  print('Completed.')
def installRMouseGeneSet(prop):
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installBMPkg.R')
  print('Install clusterProfiler (R) ...')
  assert common.runRScript(scrpt, args=['org.Mm.eg.db'], output = None)[0]
  print('Completed.')

def installEA(prop):
  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installBMPkg.R')
  print('Install clusterProfiler (R) ...')
  assert runRScript(scrpt, args=['clusterProfiler'], output = None)[0]
  print('Install pathview (R) ...')
  assert runRScript(scrpt, args=['pathview'], output = None)[0]
  print('Completed.')

# Graphics
def installRChart(prop):
  

  scrpt = os.path.join(os.environ['HYM_SCRIPT'], 'R', 'installPkg.R')

  print('Install ggplot2, (R) ...')
  assert common.runRScript(scrpt, args=['ggplot2', 'ggarchery'], output = None)[0]
  
  print('Completed.')
  

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
  os.chdir(os.environ['HYM_TEMP'])
  print('  Downloads sources ...') 
  assert common.execCmd('wget https://github.com/macs3-project/MACS/archive/refs/tags/'+ver+'.tar.gz', verbose=True)[0], 'Download error.'
  assert common.execCmd('tar -xvzf '+ver+'.tar.gz', verbose=True)[0], 'Expansion error.'
  os.chdir(os.path.join(os.environ['HYM_TEMP'],'MACS-'+ver[1:]))
  assert common.execCmd("sed -i \"s/numpy_requires = '>=1.17'/numpy_requires = '==1.22.4'/\" ./setup.py", verbose=True)[0]
  assert common.execCmd("pip install -e .", verbose=True)[0]
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('rm '+ver+'.tar.gz', verbose=True)[0]
  assert common.execCmd('rm -r MACS*', verbose=True)[0]  
  print('Completed.')
  print('>ver.', checkVerMACS(cfg))

# MEME
def checkMEME(cfg):
  return os.path.exists(os.path.join(os.environ['HYM_APP'],'meme/bin/meme'))
def checkVerMEME(cfg):
  ret = common.execCmd(os.path.join(os.environ['HYM_APP'],'meme/bin/meme -version'), showcmd=False)
  assert ret[0], 'MEME Suite is not installed.'
  ver = ret[1].strip()
  return ver
def installMEME(cfg, ver):
  print('Install MEME Suite ...')
  os.chdir(os.environ['HYM_TEMP'])
  print('  Downloads sources ...') 
  assert common.execCmd('wget https://meme-suite.org/meme/meme-software/'+ver+'/meme-'+ver+'.tar.gz', verbose=True)[0]
  assert common.execCmd('tar -xvzf ./meme-'+ver+'.tar.gz', verbose=True)[0]
  os.chdir('meme-'+ver)
  assert common.execCmd('sudo cpan force install Sys::Info',verbose=True)[0]
  assert common.execCmd('./configure --prefix='+os.path.join(os.environ['HYM_APP'], 'meme --enable-build-libxml2 --enable-build-libxslt'), verbose=True)[0]
  assert common.execCmd('make', verbose=True)[0]
  #assert common.execCmd('make test', verbose=True)[0]
  assert common.execCmd('sudo make install', verbose=True)[0]
  os.chdir(os.environ['HYM_TEMP'])
  assert common.execCmd('rm -r ./meme*', verbose=True)[0]
  os.chdir(os.environ['HYM_WS'])
  print('Completed.')
  print('>ver.', checkVerMEME(cfg))

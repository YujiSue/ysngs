import os
import subprocess
from subprocess import PIPE
import common

SOFTWARE_INFO = {
  'SRA': { 'ver' : '2.11.3' },
  'BWA': { 'ver':'0.7.17' },
  'TVC': { 'ver':'5.12.1' },
  'SAMTools': {'ver':'1.14' },
  'BCFTools': {'ver':'1.14' },
  'Picard': {'ver':'2.26.9' },
  'GATK': {'ver':'4.2.4.0' },
  'GDV' : {'ver':'1.3.0'},
  'Bowtie2': {'ver':'2.4.4' },
  'STAR': {'ver':'2.7.9a' },
  'Cuff': {'ver':'2.2.1' },
  'MACS': {'ver':'2.2.7.1' },
  'MEME':{'ver':'5.4.1' }
}
def version(name):
  return SOFTWARE_INFO[name]['ver']

def checkSRA():
  return os.path.exists(common.APPS_DIR + '/sra/sratools.' + SOFTWARE_INFO['SRA']['ver'])

def installSRA():
  print('Install SRA-Toolkit  ...')
  os.chdir(common.TEMPORAL)
  print('  Start downloading binary.') 
  proc = subprocess.run('wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/' + SOFTWARE_INFO['SRA']['ver'] + '/sratoolkit.' + SOFTWARE_INFO['SRA']['ver'] + '-ubuntu64.tar.gz', shell=True)
  proc = subprocess.run('tar xvf ./sratoolkit.' + SOFTWARE_INFO['SRA']['ver'] + '-ubuntu64.tar.gz', shell=True)
  if proc.returncode == 0:
    print('  Binary files expansion completed.') 
  else:
    print('  Binary files expansion failed.')
    return
  os.makedirs(common.APPS_DIR + '/sra',exist_ok=True)
  os.environ['PATH'] += ':' + common.APPS_DIR + '/sra'
  proc = subprocess.run('cp ./sratoolkit*/bin/* '+common.APPS_DIR+'/sra', shell=True)
  proc = subprocess.run('rm -r ./sratoolkit*', shell=True)
  proc = subprocess.run('vdb-config --interactive', shell=True)
  proc = subprocess.run('test-sra | grep "NCBI SRA Toolkit release version:" | head -n 1', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  os.chdir(common.WORK_SPACE)
  print('> ', proc.stdout.splitlines()[0])

def checkCut():
  proc = subprocess.run('which cutadapt', stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    return os.path.exists(proc.stdout.splitlines()[0])

def installCut():
  print('Install cutadapt  ...')
  proc = subprocess.run('pip install cutadapt', stdout=PIPE, stderr=PIPE, text=True, shell=True)
  print('Completed.')
  proc = subprocess.run('cutadapt --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('> ', proc.stdout.splitlines()[0])
  
def checkBWA():
  return os.path.exists(common.APPS_DIR + '/bwa')

def installBWA():
  print('Install BWA  ...')
  os.chdir(common.TEMPORAL)
  print('  Download source codes ...') 
  proc = subprocess.run('wget https://jaist.dl.sourceforge.net/project/bio-bwa/bwa-'+SOFTWARE_INFO['BWA']['ver']+'.tar.bz2', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/bwa-'+SOFTWARE_INFO['BWA']['ver']+'.tar.bz2'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('tar xvf ./bwa-'+SOFTWARE_INFO['BWA']['ver']+'.tar.bz2', shell=True)
  if proc.returncode == 1:
    print('  Failed to expand sources.')
    return
  os.chdir(common.TEMPORAL+'/bwa-'+SOFTWARE_INFO['BWA']['ver'])
  proc = subprocess.run('make -j8', shell=True)
  if proc.returncode == 0:
    proc = subprocess.run('cp bwa '+common.APPS_DIR, shell=True)
  if proc.returncode == 0:
    print('  Compile completed.') 
  else:
    print('  Compile failed.')
    return
  proc = subprocess.run('rm -r '+common.  TEMPORAL+'/bwa*', shell=True)
  os.chdir(common.WORK_SPACE)
  print('Installed.')
  proc = subprocess.run('bwa', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  lines = proc.stderr.splitlines()
  for line in lines:
    if line != '':
      print('> ', line)
    if line.startswith('Version'):
      break

def installArmadillo4TVC(dir):
  os.chdir(dir)
  proc = subprocess.run('wget http://updates.iontorrent.com/updates/software/external/armadillo-4.600.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('tar xvzf armadillo-4.600.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  os.chdir('armadillo-4.600.1')
  proc = subprocess.run("sed -i 's:^// #define ARMA_USE_LAPACK$:#define ARMA_USE_LAPACK:g' include/armadillo_bits/config.hpp", shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run("sed -i 's:^// #define ARMA_USE_BLAS$:#define ARMA_USE_BLAS:g'     include/armadillo_bits/config.hpp", shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('cmake .', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)

def installBamtools4TVC(dir):
  os.chdir(dir)
  proc = subprocess.run('wget updates.iontorrent.com/updates/software/external/bamtools-2.4.0.20150702+git15eadb925f.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('tar xvzf bamtools-2.4.0.20150702+git15eadb925f.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  os.makedirs('bamtools-2.4.0.20150702+git15eadb925f-build')
  os.chdir('bamtools-2.4.0.20150702+git15eadb925f-build')
  proc = subprocess.run('cmake ../bamtools-2.4.0.20150702+git15eadb925f -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  
def installHTSlib4TVC(dir):
  os.chdir(dir)
  proc = subprocess.run('wget --no-check-certificate https://github.com/samtools/htslib/archive/1.2.1.tar.gz -O htslib-1.2.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('tar xvzf htslib-1.2.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('ln -s htslib-1.2.1 htslib', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  os.chdir('htslib-1.2.1')
  proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)

def installSamtools4TVC(dir, dest):
  os.chdir(dir)
  proc = subprocess.run('wget --no-check-certificate https://github.com/samtools/samtools/archive/1.2.tar.gz -O samtools-1.2.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('tar xvzf samtools-1.2.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  os.chdir('samtools-1.2')
  proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('cp ./samtools '+dest, shell=True, stdout=PIPE, stderr=PIPE, text=True)
  
def checkTVC():
  return os.path.exists(common.APPS_DIR+'/TVC/bin/tvc')

def installTVC():
  print('Install Torrent Suite VariantCaller ...')
  os.chdir(common.TEMPORAL)
  print('  Downloads sources ...') 
  proc = subprocess.run('wget updates.iontorrent.com/tvc_standalone/tvc-'+SOFTWARE_INFO['TVC']['ver']+'.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+'.tar.gz'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  build = common.TEMPORAL+'/temp'
  os.makedirs(build,exist_ok=True)
  proc = subprocess.run('cp '+common.TEMPORAL+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+'.tar.gz '+build+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+'.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  install = common.APPS_DIR+'/TVC'
  os.makedirs(install,exist_ok=True)
  os.makedirs(install+'/bin',exist_ok=True)
  installArmadillo4TVC(build)
  print('  armadillo installed.') 
  installBamtools4TVC(build)
  print('  Bamtools installed.') 
  installHTSlib4TVC(build)
  print('  HTSlib installed.') 
  installSamtools4TVC(build, install+'/bin/')
  print('  Samtools(For TVC ver) installed.') 
  os.chdir(build)
  proc = subprocess.run('tar xvzf ./tvc-'+SOFTWARE_INFO['TVC']['ver']+'.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('  Source files expanded.') 
  os.makedirs(build+'/build',exist_ok=True)
  os.chdir(build+'/build')
  proc = subprocess.run('cmake '+build+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+' -DCMAKE_INSTALL_PREFIX:PATH='+install+' -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('  cmake buid.') 
  proc = subprocess.run('make -j8 install', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('  Compile completed.') 
  os.environ['PATH'] += ':'+common.APPS_DIR+'/TVC/bin'
  proc = subprocess.run('cp '+build+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+'/share/TVC/examples/example1/reference* '+common.REFERENCE_DIR, shell=True)
  proc = subprocess.run('cp '+build+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+'/share/TVC/examples/example1/test* '+common.TEST_DIR, shell=True)
  proc = subprocess.run('cp '+build+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+'/share/TVC/pluginMedia/configs/* '+common.PREFERENCE_DIR, shell=True)
  proc = subprocess.run('cp '+build+'/tvc-'+SOFTWARE_INFO['TVC']['ver']+'/share/TVC/sse/* '+common.PREFERENCE_DIR, shell=True)
  proc = subprocess.run(common.APPS_DIR+'/TVC/bin/samtools bam2fq'+common.TEST_DIR+'/test.bam > '+common.TEST_DIR+'/test.fq', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('rm -r '+build, shell=True)
  proc = subprocess.run('rm '+common.TEMPORAL+'/tvc*', shell=True)
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('tvc --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Installed.')
  print('> ', proc.stdout.splitlines()[0])

def checkSamtools():
  proc = subprocess.run('which samtools', stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    return os.path.exists(proc.stdout.splitlines()[0])

def installSamtools():
  print('Install Samtools ...')
  os.chdir(common.TEMPORAL)
  print('  Downloads sources ...') 
  proc = subprocess.run('wget https://github.com/samtools/samtools/releases/download/'+SOFTWARE_INFO['SAMTools']['ver']+'/samtools-'+SOFTWARE_INFO['SAMTools']['ver']+'.tar.bz2', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/samtools-'+SOFTWARE_INFO['SAMTools']['ver']+'.tar.bz2'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('tar xvf ./samtools-'+SOFTWARE_INFO['SAMTools']['ver']+'.tar.bz2', shell=True)
  if proc.returncode == 0:
    print('  Source files expanded.') 
  else:
    print('  Failed to expand sources')
    return
  os.chdir(common.TEMPORAL+'/samtools-'+SOFTWARE_INFO['SAMTools']['ver'])
  proc = subprocess.run('./configure', shell=True)
  if proc.returncode == 0:
    proc = subprocess.run('make -j8', shell=True)
  if proc.returncode == 0:
    print('  Compile completed.') 
  else:
    print('  Compile failed.')
    return
  proc = subprocess.run('sudo make install', shell=True)
  proc = subprocess.run('cp '+common.TEMPORAL+'/samtools-'+SOFTWARE_INFO['SAMTools']['ver']+'/examples/ex1.fa '+common.REFERENCE_DIR, shell=True)
  proc = subprocess.run('gunzip '+common.TEMPORAL+'/samtools-'+SOFTWARE_INFO['SAMTools']['ver']+'/examples/ex1.sam.gz', stdout=PIPE, stderr=PIPE, shell=True)
  proc = subprocess.run('cp '+common.TEMPORAL+'/samtools-'+SOFTWARE_INFO['SAMTools']['ver']+'/examples/ex1.sam ' +common.TEST_DIR, stdout=PIPE, stderr=PIPE, shell=True)
  proc = subprocess.run('rm -r '+common.TEMPORAL+'/samtools*', shell=True)
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('samtools view -b -T '+common.REFERENCE_DIR+'/ex1.fa -o '+common.TEST_DIR+'/ex1.bam '+common.TEST_DIR+'/ex1.sam', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('samtools fastq '+common.TEST_DIR+'/ex1.bam > '+common.TEST_DIR+'/ex1.fq', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  proc = subprocess.run('samtools --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  print('> ', proc.stdout.splitlines()[0])

def checkBCFtools():
  proc = subprocess.run('which bcftools', stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    return os.path.exists(proc.stdout.splitlines()[0])

def installBCFtools():
  print('Install BCFtools ...')
  os.chdir(common.TEMPORAL)
  print('  Downloads sources ...') 
  proc = subprocess.run('wget https://github.com/samtools/bcftools/releases/download/'+SOFTWARE_INFO['BCFTools']['ver']+'/bcftools-'+SOFTWARE_INFO['BCFTools']['ver']+'.tar.bz2', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/bcftools-'+SOFTWARE_INFO['BCFTools']['ver']+'.tar.bz2'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('tar xvf ./bcftools-'+SOFTWARE_INFO['BCFTools']['ver']+'.tar.bz2', shell=True)
  if proc.returncode == 0:
    print('  Source files expanded.') 
  else:
    print('  Failed to expand sources')
    return
  os.chdir(common.TEMPORAL+'/bcftools-'+SOFTWARE_INFO['BCFTools']['ver'])
  proc = subprocess.run('./configure', shell=True)
  if proc.returncode == 0:
    proc = subprocess.run('make -j8', shell=True)
  if proc.returncode == 0:
    print('  Compile completed.') 
  else:
    print('  Compile failed.')
    return
  proc = subprocess.run('sudo make install', shell=True)
  proc = subprocess.run('rm -r '+common.TEMPORAL+'/bcftools*', shell=True)
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('bcftools --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  print('> ', proc.stdout.splitlines()[0])

def checkPicard():
  return os.path.exists(common.APPS_DIR+'/picard.jar')

def installPicard():
  print('Install Picard ...')
  print('  Download JAR ...') 
  proc = subprocess.run('wget https://github.com/broadinstitute/picard/releases/download/'+SOFTWARE_INFO['Picard']['ver']+'/picard.jar -O '+common.APPS_DIR+'/picard.jar', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.APPS_DIR+'/picard.jar'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('java -jar '+common.APPS_DIR+'/picard.jar MarkDuplicates --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  print('> ', proc.stderr)

def checkGATK():
  return os.path.exists(common.APPS_DIR+'/gatk')

def installGATK():
  print('Install GATK ...')
  os.chdir(common.TEMPORAL)
  print('  Downloads sources ...') 
  proc = subprocess.run('wget https://github.com/broadinstitute/gatk/releases/download/'+SOFTWARE_INFO['GATK']['ver']+'/gatk-'+SOFTWARE_INFO['GATK']['ver']+'.zip', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/gatk-'+SOFTWARE_INFO['GATK']['ver']+'.zip'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('unzip -o ./gatk-'+SOFTWARE_INFO['GATK']['ver']+'.zip', shell=True)
  if proc.returncode == 1:
    print('  Failed to expand sources')
    return
  os.makedirs(common.APPS_DIR+'/gatk',exist_ok=True)
  os.environ['PATH'] += ':'+common.APPS_DIR+'/gatk'
  proc = subprocess.run('mv ./gatk-'+SOFTWARE_INFO['GATK']['ver']+'/* '+common.APPS_DIR+'/gatk', shell=True)
  proc = subprocess.run('rm -r '+common.TEMPORAL+'/gatk*', shell=True)
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('gatk --list', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  print('> ', proc.stderr.splitlines()[0])

def checkBowtie():
  proc = subprocess.run('which bowtie2', stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    return os.path.exists(proc.stdout.splitlines()[0])

def installBowtie():
  print('Install Bowtie2 ...')
  os.chdir(common.TEMPORAL)
  print('  Downloads sources ...') 
  proc = subprocess.run('wget https://jaist.dl.sourceforge.net/project/bowtie-bio/bowtie2/'+SOFTWARE_INFO['Bowtie']['ver']+'/bowtie2-'+SOFTWARE_INFO['Bowtie']['ver']+'-source.zip', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/bowtie2-'+SOFTWARE_INFO['Bowtie']['ver']+'-source.zip'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('unzip -o ./bowtie2-'+SOFTWARE_INFO['Bowtie']['ver']+'-source.zip', shell=True)
  if proc.returncode == 1:
    print('  Failed to expand sources')
    return
  os.chdir(common.TEMPORAL+'/bowtie2-'+SOFTWARE_INFO['Bowtie']['ver'])
  proc = subprocess.run('make -j8', shell=True)
  if proc.returncode == 0:
    proc = subprocess.run('sudo make install', shell=True)
  if proc.returncode == 0:
    print('  Compile completed.') 
  else:
    print('  Compile failed.')
    return
  os.chdir(common.TEMPORAL)
  proc = subprocess.run('rm -r ./bowtie2*', shell=True)
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('bowtie2 --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  print('> ', proc.stdout.splitlines()[0])

def checkSTAR():
  return os.path.exists(common.APPS_DIR+'/STAR')

def installSTAR():
  print('Install STAR ...')
  os.chdir(common.TEMPORAL)
  print('  Downloads sources ...') 
  proc = subprocess.run('wget https://github.com/alexdobin/STAR/archive/'+SOFTWARE_INFO['STAR']['ver']+'.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/'+SOFTWARE_INFO['STAR']['ver']+'.tar.gz'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('tar xvf ./'+SOFTWARE_INFO['STAR']['ver']+'.tar.gz', shell=True)
  if proc.returncode == 0:
    print('  Source files expanded.') 
  else:
    print('  Failed to expand sources')
    return
  os.chdir(common.TEMPORAL+'/STAR-'+SOFTWARE_INFO['STAR']['ver']+'/source')
  proc = subprocess.run('make -j8 STAR', shell=True)
  if proc.returncode == 0:
    print('  Compile completed.') 
  else:
    print('  Compile failed.')
    return
  proc = subprocess.run('cp STAR '+common.APPS_DIR, shell=True)
  proc = subprocess.run('rm -r '+common.TEMPORAL+'/STAR*', shell=True)
  proc = subprocess.run('rm '+common.TEMPORAL+'/'+SOFTWARE_INFO['STAR']['ver']+'*', shell=True)
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('STAR', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  lines = proc.stdout.splitlines()
  for line in lines:
    if ("version" in line) :
      print('> ', line)
      break

def checkCuff():
  return os.path.exists(common.APPS_DIR+'/cuff/cufflinks')

def installCuff():
  print('Install Cufflinks ...')
  os.chdir(common.TEMPORAL)
  print('  Download binries ...') 
  proc = subprocess.run('wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-'+SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.TEMPORAL+'/cufflinks-'+SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64.tar.gz'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('tar xvf ./cufflinks-'+SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64.tar.gz', shell=True)
  if proc.returncode == 0:
    print('  Binary files expanded.') 
  else:
    print('  Binary files expansion failed.')
    return
  os.makedirs(common.APPS_DIR+'/cuff',exist_ok=True)
  os.environ['PATH'] += ':'+common.APPS_DIR+'/cuff'
  proc = subprocess.run('mv ./cufflinks-'+SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64/* '+common.APPS_DIR+'/cuff', shell=True)
  proc = subprocess.run('rm -r ./cufflinks*', shell=True)
  proc = subprocess.run('cufflinks', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  print('> ', proc.stderr.splitlines()[0])

def checkMACS():
  proc = subprocess.run('which macs2', stdout=PIPE, stderr=PIPE, text=True, shell=True)
  if proc.returncode == 0 and proc.stdout:
    return os.path.exists(proc.stdout.splitlines()[0])

def installMACS():
  print('Install MACS2 ...')
  proc = subprocess.run('pip install macs2', stdout=PIPE, stderr=PIPE, shell=True)
  os.chdir(common.WORK_SPACE)
  proc = subprocess.run('source ~/.profile', stdout=PIPE, stderr=PIPE, shell=True)
  proc = subprocess.run('macs2 --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  print('Completed.')
  print('> ', proc.stdout.splitlines()[0])

def checkMEME():
  return os.path.exists(common.APPS_DIR+'/meme/bin/meme')

def installMEME():
  print('Install MEME Suite ...')
  os.chdir(common.TEMPORAL)
  print('  Downloads sources ...') 
  proc = subprocess.run('wget https://meme-suite.org/meme/meme-software/'+SOFTWARE_INFO['MEME']['ver']+'/meme-'+SOFTWARE_INFO['MEME']['ver']+'.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
  if proc.returncode == 0 and os.path.exists(common.common.TEMPORAL+'/meme-'+SOFTWARE_INFO['MEME']['ver']+'.tar.gz'):
    print('  Completed.') 
  else:
    print('  Failed.')
    return
  proc = subprocess.run('tar zxf ./meme-'+SOFTWARE_INFO['MEME']['ver']+'.tar.gz', shell=True)
  if proc.returncode == 0:
    print('  Source files expanded.') 
  else:
    print('  Failed to expand sources')
    return
  os.chdir('./meme-'+SOFTWARE_INFO['MEME']['ver'])
  proc = subprocess.run('./configure --prefix='+common.APPS_DIR+'/meme --enable-build-libxml2 --enable-build-libxslt', shell=True)
  proc = subprocess.run('make', shell=True)
  proc = subprocess.run('make test', shell=True)
  proc = subprocess.run('sudo make install', shell=True)
  if proc.returncode == 0:
    print('  Compile completed.') 
  else:
    print('  Compile failed.')
    return
  os.environ['PATH'] += ':'+common.APPS_DIR+'/meme/bin'
  os.chdir(common.TEMPORAL)
  proc = subprocess.run('rm -r ./meme*', shell=True)
  proc = subprocess.run('meme -version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  os.chdir(common.WORK_SPACE)
  print('Completed.')
  print('> MEME Suite', proc.stdout.splitlines()[0])

def install(exe):
  if exe == 'SRA':
    hasSRA = checkSRA()
    print('Check SRA ...', 'Installed.' if hasSRA else 'Not installed.')
    if not hasSRA:
      installSRA()
  elif exe == 'BWA':
    hasBWA = checkBWA()
    print('Check BWA ...', 'Installed.' if hasBWA else 'Not installed.')
    if not hasBWA:
      installBWA()
  elif exe == 'SAM':
    hasST = checkSamtools()
    print('Check samtools ...', 'Installed.' if hasST else 'Not installed.')
    if not hasST:
      installSamtools()
  elif exe == 'BCF':
    hasBT = checkBCFtools()
    print('Check BCFtools ...', 'Installed.' if hasBT else 'Not installed.')
    if not hasBT:
      installBCFtools()
  elif exe == 'Picard':
    hasPicard = checkPicard()
    print('Check Picard ...', 'Installed.' if hasPicard else 'Not installed.')
    if hasPicard == False:
      installPicard()
  elif exe == 'GATK':
    hasGATK = checkGATK()
    print('Check GATK ...', 'Installed.' if hasGATK else 'Not installed.')
    if hasGATK == False:
      installGATK()
  elif exe == 'Bowtie':
    hasBT = checkBowtie()
    print('Check bowtie2 ...', 'Installed.' if hasBT else 'Not installed.')
    if hasBT == False:
      installBowtie()
  elif exe == 'STAR':
    hasSTAR = checkSTAR()
    print('Check STAR ...', 'Installed.' if hasSTAR else 'Not installed.')
    if hasSTAR == False:
      installSTAR()
  elif exe == 'TVC':
    hasTVC = checkTVC()
    print('Check TVC ...', 'Installed.' if hasTVC else 'Not installed.')
    if not hasTVC:
      installTVC()
  elif exe == 'Cuff':
    hasCuff = checkCuff()
    print('Check Cufflinks ...', 'Installed.' if hasCuff else 'Not installed.')
    if hasCuff == False:
      installCuff()
  elif exe == 'MACS':
    hasMACS = checkMACS()
    print('Check MACS2 ...', 'Installed.' if hasMACS else 'Not installed.')
    if hasMACS == False:
      installMACS()
  elif exe == 'MEME':
    hasMEME = checkMEME()
    print('Check MEME ...', 'Installed.' if hasMEME else 'Not installed.')
    if hasMEME == False:
      installMEME()

import os
import subprocess
from subprocess import PIPE
import common as common

class installer :
  def __init__(self, config):
    self.cfg = config
    
  def version(self, name):
    return self.cfg.SOFTWARE_INFO[name]['ver']
    
  def checkSRA(self):
    return os.path.exists(self.cfg.APPS_DIR + '/sra/sratools.' + self.cfg.SOFTWARE_INFO['SRA']['ver'])
    
  def installSRA(self):
    print('Install SRA-Toolkit  ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Start downloading binary.') 
    proc = subprocess.run('wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/' + \
      self.cfg.SOFTWARE_INFO['SRA']['ver'] + '/sratoolkit.' + \
      self.cfg.SOFTWARE_INFO['SRA']['ver'] + '-ubuntu64.tar.gz', shell=True)
    proc = subprocess.run('tar xvf ./sratoolkit.' + self.cfg.SOFTWARE_INFO['SRA']['ver'] + '-ubuntu64.tar.gz', shell=True)
    if proc.returncode == 0:
      print('  Binary files expansion completed.') 
    else:
      print('  Binary files expansion failed.')
      return
    os.makedirs(self.cfg.APPS_DIR + '/sra',exist_ok=True)
    common.addPath(self.cfg.APPS_DIR + '/sra')
    os.system('cp ./sratoolkit*/bin/* ' + os.path.join(self.cfg.APPS_DIR,'sra'))
    os.system('rm -r ./sratoolkit*')
    proc = subprocess.run('test-sra | grep "NCBI SRA Toolkit release version:" | head -n 1', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    os.chdir(self.cfg.WORK_SPACE)
    print('> ', proc.stdout.splitlines()[0])

  def checkCut(self):
    proc = subprocess.run('which cutadapt', stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode == 0 and proc.stdout:
      return os.path.exists(proc.stdout.splitlines()[0])

  def installCut(self):
    print('Install cutadapt  ...')
    proc = subprocess.run('pip install cutadapt', stdout=PIPE, stderr=PIPE, text=True, shell=True)
    print('Completed.')
    proc = subprocess.run('cutadapt --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('> ', proc.stdout.splitlines()[0])
  
  def checkFP(self):
    return os.path.exists(os.path.join(self.cfg.APPS_DIR, 'fastp'))
  
  def installFP(self):
    print('Install fastp  ...')
    proc = subprocess.run('wget http://opengene.org/fastp/fastp -O '+os.path.join(self.cfg.APPS_DIR, 'fastp'), stdout=PIPE, stderr=PIPE, text=True, shell=True)
    proc = subprocess.run('chmod a+x '+os.path.join(self.cfg.APPS_DIR, 'fastp'), stdout=PIPE, stderr=PIPE, text=True, shell=True)
    print('Completed.')
    common.addPath(self.cfg.APPS_DIR)
    proc = subprocess.run('fastp --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('> ', proc.stderr.splitlines()[0])
  
  def checkBWA(self):
    return os.path.exists(self.cfg.APPS_DIR + '/bwa')

  def installBWA(self):
    print('Install BWA  ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Download source codes ...') 
    proc = subprocess.run('wget https://jaist.dl.sourceforge.net/project/bio-bwa/bwa-' + \
      self.cfg.SOFTWARE_INFO['BWA']['ver']+'.tar.bz2', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/bwa-' + self.cfg.SOFTWARE_INFO['BWA']['ver']+'.tar.bz2'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('tar xvf ./bwa-' + self.cfg.SOFTWARE_INFO['BWA']['ver']+'.tar.bz2', shell=True)
    if proc.returncode == 1:
      print('  Failed to expand sources.')
      return
    os.chdir(self.cfg.TEMPORAL+'/bwa-' + self.cfg.SOFTWARE_INFO['BWA']['ver'])
    proc = subprocess.run('make -j8', shell=True)
    if proc.returncode == 0:
      proc = subprocess.run('cp bwa '+self.cfg.APPS_DIR, shell=True)
    if proc.returncode == 0:
      print('  Compile completed.') 
    else:
      print('  Compile failed.')
      return
    proc = subprocess.run('rm -r '+self.cfg.TEMPORAL+'/bwa*', shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    print('Installed.')
    common.addPath(self.cfg.APPS_DIR)
    proc = subprocess.run('bwa', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    lines = proc.stderr.splitlines()
    for line in lines:
      if line != '':
        print('> ', line)
      if line.startswith('Version'):
        break

  def installArmadillo4TVC(sef, dir):
    os.chdir(dir)
    proc = subprocess.run('wget http://updates.iontorrent.com/updates/software/external/armadillo-4.600.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('tar xvzf armadillo-4.600.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    os.chdir('armadillo-4.600.1')
    proc = subprocess.run("sed -i 's:^// #define ARMA_USE_LAPACK$:#define ARMA_USE_LAPACK:g' include/armadillo_bits/config.hpp", shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run("sed -i 's:^// #define ARMA_USE_BLAS$:#define ARMA_USE_BLAS:g'     include/armadillo_bits/config.hpp", shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('cmake .', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)

  def installBamtools4TVC(self, dir):
    os.chdir(dir)
    proc = subprocess.run('wget updates.iontorrent.com/updates/software/external/bamtools-2.4.0.20150702+git15eadb925f.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('tar xvzf bamtools-2.4.0.20150702+git15eadb925f.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    os.makedirs('bamtools-2.4.0.20150702+git15eadb925f-build')
    os.chdir('bamtools-2.4.0.20150702+git15eadb925f-build')
    proc = subprocess.run('cmake ../bamtools-2.4.0.20150702+git15eadb925f -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)
  
  def installHTSlib4TVC(self, dir):
    os.chdir(dir)
    proc = subprocess.run('wget --no-check-certificate https://github.com/samtools/htslib/archive/1.2.1.tar.gz -O htslib-1.2.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('tar xvzf htslib-1.2.1.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('ln -s htslib-1.2.1 htslib', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    os.chdir('htslib-1.2.1')
    proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)

  def installSamtools4TVC(self, dir, dest):
    os.chdir(dir)
    proc = subprocess.run('wget --no-check-certificate https://github.com/samtools/samtools/archive/1.2.tar.gz -O samtools-1.2.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('tar xvzf samtools-1.2.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    os.chdir('samtools-1.2')
    proc = subprocess.run('make -j8', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('cp ./samtools '+dest, shell=True, stdout=PIPE, stderr=PIPE, text=True)
  
  def checkTVC(self):
    return os.path.exists(self.cfg.APPS_DIR+'/TVC/bin/tvc')

  def installTVC(self):
    print('Install Torrent Suite VariantCaller ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Downloads sources ...') 
    proc = subprocess.run('wget updates.iontorrent.com/tvc_standalone/tvc-' + \
      self.cfg.SOFTWARE_INFO['TVC']['ver']+'.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/tvc-' + self.cfg.SOFTWARE_INFO['TVC']['ver']+'.tar.gz'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    build = self.cfg.TEMPORAL+'/temp'
    os.makedirs(build,exist_ok=True)
    proc = subprocess.run('cp '+self.cfg.TEMPORAL+'/tvc-' + \
      self.cfg.SOFTWARE_INFO['TVC']['ver']+'.tar.gz '+build+'/tvc-' + \
        self.cfg.SOFTWARE_INFO['TVC']['ver']+'.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    install = self.cfg.APPS_DIR+'/TVC'
    os.makedirs(install,exist_ok=True)
    os.makedirs(install+'/bin',exist_ok=True)
    self.installArmadillo4TVC(build)
    print('  armadillo installed.') 
    self.installBamtools4TVC(build)
    print('  Bamtools installed.') 
    self.installHTSlib4TVC(build)
    print('  HTSlib installed.') 
    self.installSamtools4TVC(build, install+'/bin/')
    print('  Samtools(For TVC ver) installed.') 
    os.chdir(build)
    proc = subprocess.run('tar xvzf ./tvc-' + self.cfg.SOFTWARE_INFO['TVC']['ver']+'.tar.gz', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('  Source files expanded.') 
    os.makedirs(build+'/build',exist_ok=True)
    os.chdir(build+'/build')
    proc = subprocess.run('cmake '+build+'/tvc-' + \
      self.cfg.SOFTWARE_INFO['TVC']['ver'] + ' -DCMAKE_INSTALL_PREFIX:PATH=' + install + ' -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('  cmake buid.') 
    proc = subprocess.run('make -j8 install', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('  Compile completed.') 
    proc = subprocess.run('cp '+build+'/tvc-' + self.cfg.SOFTWARE_INFO['TVC']['ver']+'/share/TVC/examples/example1/reference* '+self.cfg.REFERENCE_DIR, shell=True)
    proc = subprocess.run('cp '+build+'/tvc-' + self.cfg.SOFTWARE_INFO['TVC']['ver']+'/share/TVC/examples/example1/test* '+self.cfg.TEST_DIR, shell=True)
    proc = subprocess.run('cp '+build+'/tvc-' + self.cfg.SOFTWARE_INFO['TVC']['ver']+'/share/TVC/pluginMedia/configs/* '+self.cfg.PREFERENCE_DIR, shell=True)
    proc = subprocess.run('cp '+build+'/tvc-' + self.cfg.SOFTWARE_INFO['TVC']['ver']+'/share/TVC/sse/* '+self.cfg.PREFERENCE_DIR, shell=True)
    proc = subprocess.run(self.cfg.APPS_DIR+'/TVC/bin/samtools bam2fq'+self.cfg.TEST_DIR+'/test.bam > '+self.cfg.TEST_DIR+'/test.fq', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('rm -r '+build, shell=True)
    proc = subprocess.run('rm '+self.cfg.TEMPORAL+'/tvc*', shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    common.addPath(self.cfg.APPS_DIR+'/TVC/bin')
    proc = subprocess.run('tvc --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Installed.')
    print('> ', proc.stdout.splitlines()[0])

  def checkSamtools(self):
    proc = subprocess.run('which samtools', stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode == 0 and proc.stdout:
      return os.path.exists(proc.stdout.splitlines()[0])

  def installSamtools(self):
    print('Install Samtools ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Downloads sources ...') 
    proc = subprocess.run('wget https://github.com/samtools/samtools/releases/download/' + \
      self.cfg.SOFTWARE_INFO['SAMTools']['ver']+'/samtools-' + self.cfg.SOFTWARE_INFO['SAMTools']['ver']+'.tar.bz2', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/samtools-' + self.cfg.SOFTWARE_INFO['SAMTools']['ver']+'.tar.bz2'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('tar xvf ./samtools-' + self.cfg.SOFTWARE_INFO['SAMTools']['ver']+'.tar.bz2', shell=True)
    if proc.returncode == 0:
      print('  Source files expanded.') 
    else:
      print('  Failed to expand sources')
      return
    os.chdir(self.cfg.TEMPORAL+'/samtools-' + self.cfg.SOFTWARE_INFO['SAMTools']['ver'])
    proc = subprocess.run('./configure', shell=True)
    if proc.returncode == 0:
      proc = subprocess.run('make -j8', shell=True)
    if proc.returncode == 0:
      print('  Compile completed.') 
    else:
      print('  Compile failed.')
      return
    proc = subprocess.run('sudo make install', shell=True)
    proc = subprocess.run('cp '+self.cfg.TEMPORAL+'/samtools-' + self.cfg.SOFTWARE_INFO['SAMTools']['ver']+'/examples/ex1.fa '+self.cfg.REFERENCE_DIR, shell=True)
    proc = subprocess.run('gunzip '+self.cfg.TEMPORAL+'/samtools-' + self.cfg.SOFTWARE_INFO['SAMTools']['ver']+'/examples/ex1.sam.gz', stdout=PIPE, stderr=PIPE, shell=True)
    proc = subprocess.run('cp '+self.cfg.TEMPORAL+'/samtools-' + self.cfg.SOFTWARE_INFO['SAMTools']['ver']+'/examples/ex1.sam ' +self.cfg.TEST_DIR, stdout=PIPE, stderr=PIPE, shell=True)
    proc = subprocess.run('rm -r '+self.cfg.TEMPORAL+'/samtools*', shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    proc = subprocess.run('samtools view -b -T '+self.cfg.REFERENCE_DIR+'/ex1.fa -o '+self.cfg.TEST_DIR+'/ex1.bam '+self.cfg.TEST_DIR+'/ex1.sam', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('samtools fastq '+self.cfg.TEST_DIR+'/ex1.bam > '+self.cfg.TEST_DIR+'/ex1.fq', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    proc = subprocess.run('samtools --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    print('> ', proc.stdout.splitlines()[0])

  def checkBCFtools(self):
    proc = subprocess.run('which bcftools', stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode == 0 and proc.stdout:
      return os.path.exists(proc.stdout.splitlines()[0])

  def installBCFtools(self):
    print('Install BCFtools ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Downloads sources ...') 
    proc = subprocess.run('wget https://github.com/samtools/bcftools/releases/download/' + \
      self.cfg.SOFTWARE_INFO['BCFTools']['ver']+'/bcftools-' + self.cfg.SOFTWARE_INFO['BCFTools']['ver']+'.tar.bz2', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/bcftools-' + self.cfg.SOFTWARE_INFO['BCFTools']['ver']+'.tar.bz2'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('tar xvf ./bcftools-' + self.cfg.SOFTWARE_INFO['BCFTools']['ver']+'.tar.bz2', shell=True)
    if proc.returncode == 0:
      print('  Source files expanded.') 
    else:
      print('  Failed to expand sources')
      return
    os.chdir(self.cfg.TEMPORAL+'/bcftools-' + self.cfg.SOFTWARE_INFO['BCFTools']['ver'])
    proc = subprocess.run('./configure', shell=True)
    if proc.returncode == 0:
      proc = subprocess.run('make -j8', shell=True)
    if proc.returncode == 0:
      print('  Compile completed.') 
    else:
      print('  Compile failed.')
      return
    proc = subprocess.run('sudo make install', shell=True)
    proc = subprocess.run('rm -r '+self.cfg.TEMPORAL+'/bcftools*', shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    proc = subprocess.run('bcftools --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    print('> ', proc.stdout.splitlines()[0])

  def checkPicard(self):
    return os.path.exists(self.cfg.APPS_DIR+'/picard.jar')

  def installPicard(self):
    print('Install Picard ...')
    print('  Download JAR ...') 
    proc = subprocess.run('wget https://github.com/broadinstitute/picard/releases/download/' + \
      self.cfg.SOFTWARE_INFO['Picard']['ver']+'/picard.jar -O '+self.cfg.APPS_DIR+'/picard.jar', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.APPS_DIR+'/picard.jar'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    os.chdir(self.cfg.WORK_SPACE)
    proc = subprocess.run('java -jar ' + self.cfg.APPS_DIR + '/picard.jar MarkDuplicates --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    print('> ', proc.stderr)

  def checkGATK(self):
    return os.path.exists(self.cfg.APPS_DIR + '/gatk')

  def installGATK(self):
    print('Install GATK ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Downloads sources ...') 
    proc = subprocess.run('wget https://github.com/broadinstitute/gatk/releases/download/' + \
      self.cfg.SOFTWARE_INFO['GATK']['ver']+'/gatk-' + self.cfg.SOFTWARE_INFO['GATK']['ver']+'.zip', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/gatk-' + self.cfg.SOFTWARE_INFO['GATK']['ver']+'.zip'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('unzip -o ./gatk-' + self.cfg.SOFTWARE_INFO['GATK']['ver']+'.zip', shell=True)
    if proc.returncode == 1:
      print('  Failed to expand sources')
      return
    os.makedirs(self.cfg.APPS_DIR+'/gatk',exist_ok=True)
    proc = subprocess.run('mv ./gatk-' + self.cfg.SOFTWARE_INFO['GATK']['ver']+'/* '+self.cfg.APPS_DIR+'/gatk', shell=True)
    proc = subprocess.run('rm -r '+self.cfg.TEMPORAL+'/gatk*', shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    common.addPath(self.cfg.APPS_DIR+'/gatk')
    proc = subprocess.run('gatk --list', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    print('> ', proc.stderr.splitlines()[0])

  def checkBowtie(self):
    proc = subprocess.run('which bowtie2', stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode == 0 and proc.stdout:
      return os.path.exists(proc.stdout.splitlines()[0])

  def installBowtie(self):
    print('Install Bowtie2 ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Downloads sources ...') 
    proc = subprocess.run('wget https://jaist.dl.sourceforge.net/project/bowtie-bio/bowtie2/' + \
      self.cfg.SOFTWARE_INFO['Bowtie2']['ver']+'/bowtie2-' + self.cfg.SOFTWARE_INFO['Bowtie2']['ver']+'-source.zip', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/bowtie2-' + self.cfg.SOFTWARE_INFO['Bowtie2']['ver']+'-source.zip'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('unzip -o ./bowtie2-' + self.cfg.SOFTWARE_INFO['Bowtie2']['ver']+'-source.zip', shell=True)
    if proc.returncode == 1:
      print('  Failed to expand sources')
      return
    os.chdir(self.cfg.TEMPORAL+'/bowtie2-' + self.cfg.SOFTWARE_INFO['Bowtie2']['ver'])
    proc = subprocess.run('make -j8', shell=True)
    if proc.returncode == 0:
      proc = subprocess.run('sudo make install', shell=True)
    if proc.returncode == 0:
      print('  Compile completed.') 
    else:
      print('  Compile failed.')
      return
    os.chdir(self.cfg.TEMPORAL)
    proc = subprocess.run('rm -r ./bowtie2*', shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    proc = subprocess.run('bowtie2 --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    print('> ', proc.stdout.splitlines()[0])

  def checkSTAR(self):
    return os.path.exists(self.cfg.APPS_DIR+'/STAR')

  def installSTAR(self):
    print('Install STAR ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Downloads sources ...') 
    proc = subprocess.run('wget https://github.com/alexdobin/STAR/archive/' + \
      self.cfg.SOFTWARE_INFO['STAR']['ver']+'.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/' + self.cfg.SOFTWARE_INFO['STAR']['ver']+'.tar.gz'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('tar xvf ./' + self.cfg.SOFTWARE_INFO['STAR']['ver']+'.tar.gz', shell=True)
    if proc.returncode == 0:
      print('  Source files expanded.') 
    else:
      print('  Failed to expand sources')
      return
    os.chdir(self.cfg.TEMPORAL+'/STAR-' + self.cfg.SOFTWARE_INFO['STAR']['ver']+'/source')
    proc = subprocess.run('make -j8 STAR', shell=True)
    if proc.returncode == 0:
      print('  Compile completed.') 
    else:
      print('  Compile failed.')
      return
    proc = subprocess.run('cp STAR '+self.cfg.APPS_DIR, shell=True)
    proc = subprocess.run('rm -r '+self.cfg.TEMPORAL+'/STAR*', shell=True)
    proc = subprocess.run('rm '+self.cfg.TEMPORAL+'/' + self.cfg.SOFTWARE_INFO['STAR']['ver']+'*', shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    proc = subprocess.run('STAR', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    lines = proc.stdout.splitlines()
    for line in lines:
      if ("version" in line) :
        print('> ', line)
        break
  
  def CheckHTSeq(self):
    proc = subprocess.run('which htseq-count', stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode == 0 and proc.stdout:
      return os.path.exists(proc.stdout.splitlines()[0])

  def installHTSeq(self):
    print('Install HTSeq ...')
    proc = subprocess.run('pip install HTSeq', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0:
      print('Completed.') 
    else:
      print('Failed.')
    proc = subprocess.run('htseq-count --version', stdout=PIPE, stderr=PIPE, shell=True)
    print('> ', proc.stdout.splitlines()[0])

  def checkCuff(self):
    return os.path.exists(self.cfg.APPS_DIR+'/cuff/cufflinks')

  def installCuff(self):
    print('Install Cufflinks ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Download binries ...') 
    proc = subprocess.run('wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-' + self.cfg.SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.TEMPORAL+'/cufflinks-' + self.cfg.SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64.tar.gz'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('tar xvf ./cufflinks-' + self.cfg.SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64.tar.gz', shell=True)
    if proc.returncode == 0:
      print('  Binary files expanded.') 
    else:
      print('  Binary files expansion failed.')
      return
    os.makedirs(self.cfg.APPS_DIR+'/cuff',exist_ok=True)
    os.environ['PATH'] += ':'+self.cfg.APPS_DIR+'/cuff'
    proc = subprocess.run('mv ./cufflinks-' + self.cfg.SOFTWARE_INFO['Cuff']['ver']+'.Linux_x86_64/* '+self.cfg.APPS_DIR+'/cuff', shell=True)
    proc = subprocess.run('rm -r ./cufflinks*', shell=True)
    proc = subprocess.run('cufflinks', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    print('> ', proc.stderr.splitlines()[0])

  def checkMACS(self):
    proc = subprocess.run('which macs2', stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode == 0 and proc.stdout:
      return os.path.exists(proc.stdout.splitlines()[0])

  def installMACS(self):
    print('Install MACS2 ...')
    proc = subprocess.run('pip install macs2', stdout=PIPE, stderr=PIPE, shell=True)
    os.chdir(self.cfg.WORK_SPACE)
    proc = subprocess.run('source ~/.profile', stdout=PIPE, stderr=PIPE, shell=True)
    proc = subprocess.run('macs2 --version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('Completed.')
    print('> ', proc.stdout.splitlines()[0])

  def checkMEME(self):
    return os.path.exists(self.cfg.APPS_DIR+'/meme/bin/meme')

  def installMEME(self):
    print('Install MEME Suite ...')
    os.chdir(self.cfg.TEMPORAL)
    print('  Downloads sources ...') 
    proc = subprocess.run('wget https://meme-suite.org/meme/meme-software/' + \
      self.cfg.SOFTWARE_INFO['MEME']['ver']+'/meme-' + self.cfg.SOFTWARE_INFO['MEME']['ver'] + '.tar.gz', stdout=PIPE, stderr=PIPE, shell=True)
    if proc.returncode == 0 and os.path.exists(self.cfg.self.cfg.TEMPORAL+'/meme-' + self.cfg.SOFTWARE_INFO['MEME']['ver']+'.tar.gz'):
      print('  Completed.') 
    else:
      print('  Failed.')
      return
    proc = subprocess.run('tar zxf ./meme-' + self.cfg.SOFTWARE_INFO['MEME']['ver']+'.tar.gz', shell=True)
    if proc.returncode == 0:
      print('  Source files expanded.') 
    else:
      print('  Failed to expand sources')
      return
    os.chdir('./meme-' + self.cfg.SOFTWARE_INFO['MEME']['ver'])
    proc = subprocess.run('./configure --prefix='+self.cfg.APPS_DIR+'/meme --enable-build-libxml2 --enable-build-libxslt', shell=True)
    proc = subprocess.run('make', shell=True)
    proc = subprocess.run('make test', shell=True)
    proc = subprocess.run('sudo make install', shell=True)
    if proc.returncode == 0:
      print('  Compile completed.') 
    else:
      print('  Compile failed.')
      return
    os.environ['PATH'] += ':'+self.cfg.APPS_DIR+'/meme/bin'
    os.chdir(self.cfg.TEMPORAL)
    proc = subprocess.run('rm -r ./meme*', shell=True)
    proc = subprocess.run('meme -version', shell=True, stdout=PIPE, stderr=PIPE, text=True)
    os.chdir(self.cfg.WORK_SPACE)
    print('Completed.')
    print('> MEME Suite', proc.stdout.splitlines()[0])

  def install(self, exe):
    if exe == 'SRA':
      hasSRA = self.checkSRA()
      print('Check SRA ...', 'Installed.' if hasSRA else 'Not installed.')
      if not hasSRA:
        self.installSRA()
    elif exe == 'Cutter':
      hasCut = self.checkCut()
      print('Check CutAdapt ...', 'Installed.' if hasCut else 'Not installed.')
      if not hasCut:
        self.installCut()
    elif exe == 'FP':
      hasFP = self.checkFP()
      print('Check fastp ...', 'Installed.' if hasFP else 'Not installed.')
      if not hasFP:
        self.installFP()
    elif exe == 'BWA':
      hasBWA = self.checkBWA()
      print('Check BWA ...', 'Installed.' if hasBWA else 'Not installed.')
      if not hasBWA:
        self.installBWA()
    elif exe == 'SAM':
      hasST = self.checkSamtools()
      print('Check samtools ...', 'Installed.' if hasST else 'Not installed.')
      if not hasST:
        self.installSamtools()
    elif exe == 'BCF':
      hasBT = self.checkBCFtools()
      print('Check BCFtools ...', 'Installed.' if hasBT else 'Not installed.')
      if not hasBT:
        self.installBCFtools()
    elif exe == 'Picard':
      hasPicard = self.checkPicard()
      print('Check Picard ...', 'Installed.' if hasPicard else 'Not installed.')
      if not hasPicard:
        self.installPicard()
    elif exe == 'GATK':
      hasGATK = self.checkGATK()
      print('Check GATK ...', 'Installed.' if hasGATK else 'Not installed.')
      if not hasGATK:
        self.installGATK()
    elif exe == 'Bowtie':
      hasBwT = self.checkBowtie()
      print('Check bowtie2 ...', 'Installed.' if hasBwT else 'Not installed.')
      if not hasBwT:
        self.installBowtie()
    elif exe == 'STAR':
      hasSTAR = self.checkSTAR()
      print('Check STAR ...', 'Installed.' if hasSTAR else 'Not installed.')
      if not hasSTAR:
        self.installSTAR()
    elif exe == 'TVC':
      hasTVC = self.checkTVC()
      print('Check TVC ...', 'Installed.' if hasTVC else 'Not installed.')
      if not hasTVC:
        self.installTVC()
    elif exe == 'HTS':
      hasHTS = self.checkCuff()
      print('Check HTSeq ...', 'Installed.' if hasHTS else 'Not installed.')
      if not hasHTS:
         self.installHTSeq()
    elif exe == 'Cuff':
      hasCuff = self.checkCuff()
      print('Check Cufflinks ...', 'Installed.' if hasCuff else 'Not installed.')
      if not hasCuff:
         self.installCuff()
    elif exe == 'MACS':
      hasMACS = self.checkMACS()
      print('Check MACS2 ...', 'Installed.' if hasMACS else 'Not installed.')
      if not hasMACS:
        self.installMACS()
    elif exe == 'MEME':
      hasMEME = self.checkMEME()
      print('Check MEME ...', 'Installed.' if hasMEME else 'Not installed.')
      if not hasMEME:
        self.installMEME()

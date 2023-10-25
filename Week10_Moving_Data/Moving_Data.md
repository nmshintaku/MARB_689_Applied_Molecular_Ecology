# **How to transfer files to and from the TAMU HPRC**

Today, we will prepare for you semester project by downloading your chosen datasets to the server. We will also practice how to download data to your personal computer, whether they be tables or figures, for your final report. Additional details on how to transfer files on the TAMU HPRC can be found here: https://hprc.tamu.edu/wiki/HPRC:File_Transfers.

## Getting data on the server
### Moving files from your computer on and off the server
Using a file transfer software with SCP or SFTP protocols. The software package options offer an interactive interface like your computer's file storage system that makes the file transfer process more intuitive. The files are securely transferred between the local host and remote host. This option is not the fastest because it is working directly on the login node to transfer the files.

In all options below you will want to set your remote host to either grace.hprc.tamu.edu or the transfer node grace-dtn1.hprc.tamu.edu. The first option will time out after 60 mins, meaning if your transferring larger files you should use the second transfer node option.

Windows Options:
A. [WinSCP](https://winscp.net/eng/index.php)
B. [FileZilla](https://filezilla-project.org/)

Mac OS Options:
A. [CyberDuck](https://cyberduck.io/)
B. [FileZilla](https://filezilla-project.org/)
C. scp/sftp commands directly from the terminal

Copy files from your computer to the server
```
scp local_filename username@grace.hprc.tamu.edu:/scratch/group/kitchen-group/class_working_directories/
```
Copy files from the server to your computer
```
scp username@grace.hprc.tamu.edu:/scratch/group/kitchen-group/class_working_directories/remote_file .
```
### Moving files to remote storage servers
Because the HPRC is not a permanent storage solution, you will need to periodically move data from the server to a cloud storage sites like Google Drive, Dropbox, BOX, Amazon AWS, etc. This can be done using *rclone* command, available on grace and terra transfer nodes servers(grace-dtn1.hprc.tamu.edu and terra-ftn.hprc.tamu.edu).

You will need to configure your destination (i.e., Google Drive). Follow instructions here: https://hprc.tamu.edu/wiki/SW:rclone

Once your have configured the drive, you can use rclone as follows.
To copy:
```
rclone copy $HOME mygdrive:home-grace
```
To sync:
```
rclone sync $HOME mygdrive:home-grace
```

### Downloading files from an online repository to the server
Commands to know:
GNU *wget* command is used to download files from the web. It uses HTTP, HTTPS and FTP protocols. The file will download to your current working directory. More info [here](https://linuxize.com/post/wget-command-examples/).
```
wget [options] [url]
```

Another command line tool, *curl*, can be used to transfer data from the web.
```
curl [options] [url]
```

#### NCBI SRA database
The proper way to download the data from SRA is to use a command line tools *prefetch* and *fasterq-dump* (https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/).

Here is an example of how that command would work.
```
#!/bin/bash
#SBATCH --time=08:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=10G   # memory per CPU core
#SBATCH --job-name "fd"
#SBATCH --output "fd.log"

# load modules
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load SRA-Toolkit/3.0.3

# enable proxy to connect to the internet
module load WebProxy

# loop to download each SRA accession in SRA.list file
for SAMPLE in $(cat SRA.list); do
        prefetch --output-directory ./ $SAMPLE && \
        fasterq-dump --outdir ./ --split-files -e 8 $SAMPLE
        pigz -p 8 $SAMPLE_1.fastq $SAMPLE_2.fastq
done

```
The batch script above depends on the file `SRA.list`. This is a simple text file where each row is an SRA accession.
```
SRR25340716
SRR25655135
SRR25338495
```

Another workaround is to use https://sra-explorer.info/. You can search for your SRA accession or project ID that will generate the *curl* commands you can run on the command line prompt to directly download your data into that directory.

#### DRYAD of Zendo
For these repositories it is easiest to use *wget* or *curl* commands.

#### GitHub
Here you can download the entire repository using the command:

```
git clone [url]
```

### Transferring files from server to server
The best option to move data between servers is to use Globus Connect. This file transfer platform can transfer larger amounts of data (100+ GB) or sync directories from the remote host to the endpoint (another server or your local computer). Use endpoints for TAMU HPRC include: "TAMU terra-ftn" for Terra cluster, and "TAMU grace-dtn" for Grace cluster.

Additional information can be found at https://hprc.tamu.edu/wiki/SW:GlobusConnect.

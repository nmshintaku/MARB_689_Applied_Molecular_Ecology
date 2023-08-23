# **The Basics of Command Line**
August 23, 2023

### Let's Login

Either in the terminal (Mac) or GitBash interface (PC), please type in:

```
ssh <NETID>@grace.hprc.tamu.edu
```
If you are using Putty, you may have already saved your configuration. In that case, click on the "Open" button. *_see prior HPC Intro document for guidance_

Input your TAMU password (hidden from view) and press Enter. A Duo push notification window will appear.

### You're in!!!

### So what's the difference between terminal and command line anyway?

Generally, you'll find these two terms used interchangeably. Technically, a terminal is a software that starts and connects to a shell. A shell is your session and session environment (where things like the prompt and shortcuts might be customized). The command line is the literal line where you enter commands and the cursor blinks.

### Now what?
Let's figure out where you are. Type in `pwd` at the command prompt. What does it say?

Now let's navigate to our shared group space:
```
cd /scratch/group/kitchen-group
```
_cd= change directory_

Next, we will see what is here. Type `ll` at the command prompt. What do you see?

What do the colors signify?
Blue =
Green =
White =
Red =
Cyan =

Now type `ls`. How is this different than `ll`?

Now type `ls -ltr`

Navigate into the directory `class_working_directories`. What's in there?

Let's create a subdirectory for your class results here. Use your **NETID** as the subdirectory name.

```
mkdir kitchens
```
_mkdir= make directory_

### **VERY IMPORTANT**
Permissions: We will want to change permission on your new folder to allow course instructors (i.e., Sheila) to assist with troubleshooting code in the future.:
```
chmod g+rw kitchens
```
_chmod=change mode, change permissions on files and directories_
_g= group, + = add, -= remove, r= read access, w= write access, x= execute_

File permissions are core to the security model used by Linux systems. They determine who can access files and directories on a system and how.

What permissions does the user have on the `00_README.txt` file in the group directory?

What permissions does the group have on the `fix_permission.kitchen-group.sh` script?

What permissions does the universe have on the `fix_permission.kitchen-group.sh` script?

How would we remove permissions for the `group` to `write` to our new directories?

### Let's run a script

Let's move into your new directory you've made.

In that directory, we are going to copy over a script from several directories away (`01_theBasics`).  
```
cp ../../01_theBasics/test.sh .
```
_cp= copy fileA to this location (.)_
_. = this directory_
_../= one directory away_

**Challenge** Check what is in your directory. What permissions do you, the user, have on the item(s) in this directory?

Open the copied script using GNU `nano` text editor:

```
nano test.sh
```

#### Components of the script
```
#!/bin/bash
#SBATCH --time=0:10:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1M   # memory per CPU core
#SBATCH --output "test_1"
```
1. The Shebang command (`#!/bin/bash`) tells the shell (which interprets the UNIX commands) to interpret and run the Slurm script using the bash (Bourne-again shell) shell. This line should always be added at the very top of your SBATCH/Slurm script.

2. In the next section, the amount of resources required for the job to run on the compute nodes are specified. This informs Slurm about the name of the job, output filename, amount of RAM, Nos. of CPUs, nodes, tasks, time, and other parameters to be used for processing the job.

These SBATCH commands are also know as SBATCH directives and must be preceded with a pound sign and should be in an uppercase format as shown above

|Specification	|Option	|
|-------------|-----|
|Wall Clock Limit	|--time=[hh:mm:ss]|
|Job Name	|--job-name=[SomeText]	|
|Total Task/Core Count	|--ntasks=[#]	|
|Tasks per Node I	|--ntasks-per-node=#	|
|Memory Per Node	|--mem=value[K\|M\|G\|T]	|
|Combined stdout/stderr	|--output=[OutputName].%j	|

3. Load any dependencies, if needed (not in this example).

4. The jobs to run:

```
echo "Time: $(date)"
start=$(date +%s)

# Change the information below
USER=kitchens
DIRNAME=bananas

mkdir /scratch/group/kitchen-group/class_working_directories/${USER}/${DIRNAME}

echo "Test Done !!!!!!!"

echo "Time: $(date)"
end=$(date +%s)
runtime=$((end-start))

echo "Runtime was $runtime"
```

In the script, we will change the variables for `USER` and `DIRNAME`.

Save your changes (ctrl + o, enter) and close (ctrl + x) the nano window.

Now its time to execute the script. This is as simple as"

```
sbatch test.sh
```

This script will finish before we can even check its run status. Let's try anyways. Type in `squeue -u <NETID>`.

Recap on important commands to submit jobs on the cluster.

|Commands	|Syntax	|Description|
|--|--|
|sbatch | sbatch <job-id>|Submit a batch script to Slurm for processing.|
|squeue | squeue -u | Show information about your job(s) in the queue. The command when run without the -u flag, shows a list of your job(s) and all other jobs in the queue.|
|skill/scancel|scancel <job-id>|End or cancel a queued job.|

#### Additional useful information on job submission on the Grace Cluster
Documentation on all things relevant to using the Grace HPRC are provided at https://hprc.tamu.edu/wiki/Grace:Batch.

### Let's practice a couple more basic commands
Navigate to the `01_theBasics` directory. What files do you see in there?

Let's take a peek at one of the zipped sequence file from Hawaiian Red Shrimp _Halocaridina rubra_, Hal_SL18084_R1_2mil.fastq.gz.
<p align="center"> <img src="./images/halocaridina-rubra.jpg" width="250" />

```
zcat Hal_SL18084_R1_2mil.fastq.gz | head -n 15
```
_zcat= look at file without uncompressing it_
_head = print N top number of lines, default is 10 lines_
_tail = like head, but the bottom N lines_
_| = called a pipe, allows you to string together two commands like zcat and head in this example_

**DO NOT RUN:** To completely uncompress the file we could run `gunzip filename` and to recompress it `gzip filename`. These are not the only compression tool options. You may encounter .zip, tar.gz or .bz files in the future. That sounds like a future problem to solve, right?

`zcat` is a version of `cat`, or concatenate. This command will read the data from the file and output it to screen. It can also be used to join two or more files from top to bottom. For example, `cat fileA fileB > fileAB`.

The `>` will write out whereas the `<` will read in to the command you are writing. See example above of how we combined two file into a new file called `fileAB`.

It can be very useful to quickly screen a file for a string or pattern. To do this we will use `grep`.

Count the number (-c) of lines "group" appears in `test.sh`?
```
grep -c "group" test.sh
```

Count the number (-c) of times "group" appears in `test.sh`?
```
grep -o "group" test.sh | wc -l
```
_wc -l = word count provides number of lines, word count, byte and characters count, adding `-l` only prints number of lines_

Extract one line prior to the first appearance of "USER".
```
grep -B 1 "USER" test.sh
```
_B= number of lines before string_
_A = number of lines after string_

### DANGER!!! Avoid these if you can.
Occasionally we want to move (`mv`) a file from one directory to another or we may want to remove (`rm`) it altogether. These commands are dangerous because they cannot be undone. Once a file is removed, it cannot be recovered. Or, if you move a file into your directory with the same name, the prior file will be overwritten. So, in most cases in this class think EXTRA EXTRA HARD about using these commands. Is it worth it?

### Creating short-cuts (i.e., aliases)
1. Navigate back home:
```
cd ~
```
OR

```
cd
```
In your home directory, type:
```
nano .bashrc
```

2. Under 'User specific aliases and functions', add:
```
alias class='cd /scratch/group/kitchen-group/class_working_directories/<NETID>'
```

3. Save your changes (ctrl + o, enter) and close (ctrl + x) the nano window.

** NOTE: This will not become active until you log out and log back in **

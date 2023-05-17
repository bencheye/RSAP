# RSAP: a RNASeq analysis pipeline

## Quick start

### preparation for the running environment

1. install docker engine

   Download and install the docker engine on the basis of [page](https://docs.docker.com/desktop/). Then create an account and sign in the account in the docker desktop.

2. download the RSAP docker image

   Download the RSAP docker image with the following code block in terminal.

   ```shell
   docker pull bencheye/rsap:tagname
   ```

3. run the RSAP docker images 

   Building the RSAP container with the RSAP image.

   ```shell
   docker run -it --name containerName -v hostPath:containerPath imageID
   # containerName refers to container name that you can customize to set.
   # hostPath refers to the host directory to save analysis project interactive 
   # with container.
   # containerPath refers to container directory interactive with host.
   ```

4. load the RSAP container 

   Loading the RSAP container  to enter the RASP running systems.

   ```shell
   # 1. open a new bash terminal (recommend)
   docker exec -it containerName /bin/bash
   # 2. attach a opening container
   docker attach containerName
   ```

5. other common docker command

   ```shell
   docker images # View all downloaded images
   docker rmi imageID/imageName # remove image
   docker ps # Show containers that are currently running
   docker ps -a # Show all the containers including stopped container
   docker rm containerName # remove container
   docker start containerName # start container
   docker stop containerName # stop container
   ```

### Run RASP 

1. Initialization settings (first running)

   After entering the container system, you can run `initiationSetup` commands anywhere to initialize RSAP running settings.`initiationSetup` must be run when first enter the RSAP container. And  when you enter the container again, it is unnecessary to run `initiationSetup` after first run it.

   ```shell
   initiationSetup workspacePath
   # workspacePath refers to container dircetory to save RSAP analysis projects.
   # it is best to be consistent with containerPath that is convenient to exchange 
   # data with host.
   # workspacePath must be an absolute path.
   ```

2. Load running environment 

   Activating the running environment required for RSAP operation.

   ```shell
   conda activate RSAP
   ```

3. Generate an analysis project

   Then, to generate an analysis project with `generateProject`. This command has two parameters. The first parameter is project name to specify a project name. Another is used to specify the sequence type, pair (paired-end) and single (single-end). After running `generateProject`,  a directory `{projectName}_RNAseq_analysis_{todayDate}_{nowTime}` will be created. And in this directory, it will generate a `metadata.xls` that is a sample mapping information between samples and groups, `project_run_config.yaml` including all the parameters needed for RSAP running and `RSAP_running.py` the RSAP running script. Meanwhile, it also generates four folders, `Data` to place original fastq sequence data, `Result` to place all the analysis results, `Log` to place RSAP running log file and `Report` to place the final report.

   ```shell
   generateProject projectName pair/single
   ```

4. Fill in personalized project information

   Firstly, please fill in the `metadata.xls` according to the original analysis data (fastq). The first column refers to `Sample` name, and the second is `Group` name. Then, please modify the configuration file, `project_run_config.yaml`, based on the operation requirements and sample characteristics.

   **Notice:** for single-end data, the sample name must be `{sample}.{Postfix}`, e.g. `A1.fastq.gz`. For paired-end data, the sample name must be `{sample}_{PairedSep}.{Postfix}`, e.g. `A1_R1.fastq.gz`. `Sample` and `Group` can be generated automatically according to the `metadata.xls` without filling.

5. Run RSAP

   Finally, run the RSAP by `python RSAP_running.py`.

   


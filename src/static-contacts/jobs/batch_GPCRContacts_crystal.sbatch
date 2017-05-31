#!/bin/bash

#Sbatch for: Computing all GPCRContacts for class A crystal structures
#################
#set a job name
#SBATCH --job-name=GPCRContacts_crystal
#SBATCH --mail-user=akma327@stanford.edu --mail-type=ALL
#################
# A file for job output, you can check job progress
#SBATCH --output=/scratch/PI/rondror/akma327/GPCRContacts/src/static-contacts/jobs/GPCRContacts_crystal_PROGRESS.out
#################
# A file for errors from the job
#SBATCH --error=/scratch/PI/rondror/akma327/GPCRContacts/src/static-contacts/jobs/GPCRContacts_crystal_ERROR.out
#################
#SBATCH --time=5:00:00
#################
#SBATCH --partition=rondror
#SBATCH --qos=rondror
#################
#number of nodes you are requesting
#SBATCH --tasks=4
#SBATCH --ntasks-per-socket=6
#SBATCH --mem=60000
#################

cd /scratch/PI/rondror/akma327/GPCRContacts/src/static-contacts
python batch_static_contacts_classA_gpcrs.py



wait

#!/bin/bash
#Sbatch for: Generate evobundle-simulation flareplots
#################
#set a job name
#SBATCH --job-name=gen_flareplot_dataset_batch_lwb
#SBATCH --mail-user=akma327@stanford.edu --mail-type=ALL
#################
# A file for job output, you can check job progress
#SBATCH --output=gen_flareplot_dataset_batch_lwb_PROGRESS.out
#################
# A file for errors from the job
#SBATCH --error=gen_flareplot_dataset_batch_lwb_ERROR.out
#################
#SBATCH --time=1:00:00
#################
#SBATCH --partition=rondror
#SBATCH --qos=rondror
#################
#number of nodes you are requesting
#SBATCH --tasks=4
#SBATCH --ntasks-per-socket=6
#SBATCH --mem=10000
#################


EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/B2AR-inactive-carazalol-cell2013b/condition-12/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="2RH1"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/B2AR-inactive-carazalol-cell2013b-cond12-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/B2AR-active-Nb-BI-pnas2011b/condition-Q/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3P0G"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/B2AR-active-Nb-BI-pnas2011b-condQ-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/B2AR-active-Nb-BI-pnas2011b/condition-Q/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3P0G"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/B2AR-active-Nb-BI-pnas2011b-condQ-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/B2AR-active-Gs-BI-science2015/condition-22/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3SN6"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/B2AR-active-Gs-BI-science2015-cond22-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/B2AR-active-Gs-BI-science2015/condition-23/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3SN6"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/B2AR-active-Gs-BI-science2015-cond23-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/B2AR-active-Gs-BI-science2015/condition-24/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3SN6"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/B2AR-active-Gs-BI-science2015-cond24-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/B2AR-active-Gs-BI-science2015/condition-65/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3SN6"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/B2AR-active-Gs-BI-science2015-cond65-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/D4R-inactive-nemanopride-unpublished/condition-antagonist-bound/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="D4DR"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/D4R-inactive-nemanopride-unpublished-cond-antagonist-bound-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/D4R-inactive-nemanopride-unpublished/condition-antagonist-bound/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="D4DR"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/D4R-inactive-nemanopride-unpublished-cond-antagonist-bound-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/D4R-inactive-nemanopride-unpublished/condition-antagonist-bound/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="D4DR"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/D4R-inactive-nemanopride-unpublished-cond-antagonist-bound-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/DOR-inactive-naltrindole-unpublished/condition-naltrindole-bound/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4N6H"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/DOR-inactive-naltrindole-unpublished-condition-naltrindole-bound-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/DOR-inactive-naltrindole-unpublished/condition-naltrindole-bound/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4N6H"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/DOR-inactive-naltrindole-unpublished-condition-naltrindole-bound-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/DOR-inactive-naltrindole-unpublished/condition-naltrindole-bound/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4N6H"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/DOR-inactive-naltrindole-unpublished-condition-naltrindole-bound-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-inactive-QNB-nature2013/condition-AE/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3UON"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-inactive-QNB-nature2013-condAE-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-inactive-QNB-nature2012/condition-A/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="3UON"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-inactive-QNB-nature2012-condA-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-inactive-bFNA-nature2015/condition-bFNA-bound/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4DKL"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-inactive-bFNA-nature2015-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-inactive-bFNA-nature2015/condition-bFNA-bound/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4DKL"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-inactive-bFNA-nature2015-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-inactive-bFNA-nature2015/condition-bFNA-bound/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4DKL"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-inactive-bFNA-nature2015-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-active-noNb-BU72-nature2015/condition-nb-interface-restrained-BU72/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="5C1M"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-active-noNb-BU72-nature2015-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-active-noNb-BU72-nature2015/condition-nb-interface-restrained-BU72/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="5C1M"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-active-noNb-BU72-nature2015-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-active-noNb-BU72-nature2015/condition-nb-interface-restrained-BU72/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="5C1M"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-active-noNb-BU72-nature2015-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-active-Nb-unpublished/condition-prot-D114-D164-nb-interface-restrained/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="5C1M"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-active-Nb-unpublished-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-active-Nb-unpublished/condition-prot-D114-D164-nb-interface-restrained/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="5C1M"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-active-Nb-unpublished-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-active-Nb-unpublished/condition-prot-D114-D164-nb-interface-restrained/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="5C1M"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/MOR-active-Nb-unpublished-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-active-restrained-unpublished/condition-prot-asp2x50-3x49/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4MQS"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-active-restrained-unpublished-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-active-restrained-unpublished/condition-prot-asp2x50-3x49/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4MQS"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-active-restrained-unpublished-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-active-restrained-unpublished/condition-prot-asp2x50-3x49/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4MQS"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-active-restrained-unpublished-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-active-unrestrained-unpublished/condition-prot-asp2x50-3x49/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4MQS"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-active-unrestrained-unpublished-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-active-unrestrained-unpublished/condition-prot-asp2x50-3x49/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4MQS"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-active-unrestrained-unpublished-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/M2R-active-unrestrained-unpublished/condition-prot-asp2x50-3x49/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="4MQS"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/M2R-active-unrestrained-unpublished-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/Rhodopsin-active-unrestrained-unpublished/2x72-E113Q-M257Y-constitutively-active/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP1"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/RHO-active-unrestrained-unpublished-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/Rhodopsin-active-unrestrained-unpublished/2x72-E113Q-M257Y-constitutively-active/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP1"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/RHO-active-unrestrained-unpublished-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/Rhodopsin-active-unrestrained-unpublished/2x72-E113Q-M257Y-constitutively-active/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP1"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/RHO-active-unrestrained-unpublished-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/US28-active-unrestrained-unpublished/fractalkine-bound/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP2"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/US28-active-unrestrained-unpublished-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/US28-active-unrestrained-unpublished/fractalkine-bound/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP2"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/US28-active-unrestrained-unpublished-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/US28-active-unrestrained-unpublished/fractalkine-bound/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP2"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/US28-active-unrestrained-unpublished-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/US28-active-restrained-unpublished/fractalkine-bound/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP2"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/US28-active-restrained-unpublished-rep_1.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/US28-active-restrained-unpublished/fractalkine-bound/rep_2/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP2"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/US28-active-restrained-unpublished-rep_2.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &

EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/US28-active-restrained-unpublished/fractalkine-bound/rep_3/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="TMP2"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/lwb/US28-active-restrained-unpublished-rep_3.json"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -lwb &



wait

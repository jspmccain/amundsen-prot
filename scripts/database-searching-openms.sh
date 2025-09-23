#!/bin/bash
# Database searching and fdr application

database_fasta_file=$1
mzml_folder=$2

DIR=$mzml_folder
for FILE in "$DIR"*.mzML
do
	echo "Processing $FILE file..."
        temp_string=${FILE/.mzML/}

# formatting the input names so that they can properly feed into the database search
        db_string=${database_fasta_file/.fasta/}
        db_string_adjusted=$db_string'.fasta'
        db_string_revcat=$db_string'.revCat.fasta'

        echo $db_string
        echo $db_string_adjusted
        echo $db_string_revcat

# running the database search
	MSGFPlusAdapter -in $FILE -executable /software/MSGFPlus/MSGFPlus_Oct2018.jar -database $db_string_adjusted -out $temp_string'.idXML' -add_decoys -threads 30 -java_memory 50000 -fixed_modifications 'Carbamidomethyl (C)' -variable_modifications 'Oxidation (M)' 'Gln->pyro-Glu (N-term Q)' 'Deamidated (N)' 'Deamidated (Q)'
        PeptideIndexer -in $temp_string'.idXML' -fasta $db_string_revcat -out $temp_string'_PI.idXML' -threads 6 -decoy_string 'XXX_' -enzyme:specificity 'none'
        FalseDiscoveryRate -in $temp_string'_PI.idXML' -out $temp_string'_FDR.idXML' -PSM 'true' -FDR:PSM 0.01 -threads 4

done



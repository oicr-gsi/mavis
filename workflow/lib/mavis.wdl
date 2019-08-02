version 1.0

workflow MavisWorkflow {
 File inputBAM
 File inputBAMindex
 File STARFusion
 String projectID
 String outputDIR
 String outputCONFIG

 call configureMavis { input: outputDIR = outputDIR, outputCONFIG = outputCONFIG, projectID = projectID, STARFusion = STARFusion, inputBAM = inputBAM, inputBAMidx = inputBAMindex }
 call setupMavis { input: inputAnchor = configureMavis.output_message, outputDIR = outputDIR, outputCONFIGfile = configureMavis.mavis_config, STARFusion = STARFusion }
 call launchMavis { input: inputAnchor = setupMavis.output_message, outputDIR = outputDIR, outputCONFIG = outputCONFIG }
 call getJobID { input: inputStrings = launchMavis.output_message }
 call finishMavis { input: lastJobId = getJobID.job_id, outputDIR = outputDIR}

}

# ===================================
#            CONFIGURE
# ===================================
task configureMavis {

File   inputBAM
File   inputBAMidx
File   STARFusion
String outputDIR
String outputCONFIG
String projectID

command {
 set -euo pipefail
 module load python-gsi/3.6.4
 mavis config --write "${outputDIR}${outputCONFIG}" \
              --assign ${projectID} starfusion starfusion \
              --library ${projectID} transcriptome diseased True ${inputBAM} \
              --convert starfusion ${STARFusion} starfusion
 echo "MAVIS configured"
}

runtime {
  memory: "10 GB"
}

output {
  String output_message = read_string(stdout())
  File mavis_config = "${outputDIR}${outputCONFIG}"
}
}


# ====================================
#          SETUP
# ====================================
task setupMavis {

String outputDIR
String outputCONFIGfile
File   STARFusion
String inputAnchor

command {
  module load python-gsi/3.6.4
  set -euo pipefail
  mavis setup ${outputCONFIGfile} -o ${outputDIR}
  echo "MAVIS setup complete"
}

runtime {
  memory: "10 GB"
}

output {
 String output_message = read_string(stdout())
}

}

# ====================================
#          LAUNCH
# ====================================
task launchMavis {

String outputDIR
String outputCONFIG
String inputAnchor


command {
 echo "Launching Mavis"
 set -euo pipefail
 module load python-gsi/3.6.4
 mavis schedule -o ${outputDIR} --submit
}

runtime {
  memory: "12 GB"
  continueOnReturnCode: [0, 1]
}

output {
 Array[String] output_message = read_lines(stderr())
}

}

# ===========================================
#              GET ID
# ===========================================
task getJobID {
 Array[String] inputStrings

 command {
     cat ${write_lines(inputStrings)} | grep SUBMITTED | tail -n 1 | sed s/.*\(// | sed s/\).*//    
 }
 output {
     Int job_id = read_int(stdout())
 }
}

# ===========================================
#              FINISH
# ===========================================
task finishMavis {

Int lastJobId
Int? sleepInterval = 10
String outputDIR

command {
 while qstat -j ${lastJobId}; do
  sleep ${sleepInterval}
 done
 if [ -f ${outputDIR}summary/MAVIS-${lastJobId}.COMPLETE ]
    then
     exit 0;
 fi
 exit 1;
}

}

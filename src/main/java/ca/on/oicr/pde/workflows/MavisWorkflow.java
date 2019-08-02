package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MavisWorkflow extends OicrWorkflow {

    private final Logger logger = LoggerFactory.getLogger(MavisWorkflow.class);

    private String javaPath;
    private String seqwareDistributionJarPath;
    private String launchAndWaitScript;
    private String copyOutputScript;
    private String provisionOutScript;
    private String cromwellJarPath;
    private String cromwellHost;
    private String wdlWorkflow;
    private String wdlInputs;
    private String wdlOptions;
    private String wdlDepsZip;
    private String outputFilesDefinition;
    private boolean manualOutput = false;

    private String outputDir = "output/";
    private String tmpDir = "tmp/";

    private void init() {
        javaPath = getRequiredProperty("java_path");
        seqwareDistributionJarPath = getRequiredProperty("seqware_distribution_jar_path");
        launchAndWaitScript = getRequiredProperty("launch_and_wait_script");
        copyOutputScript = getRequiredProperty("copy_output_script");
        provisionOutScript = getRequiredProperty("provision_out_script");

        cromwellJarPath = getRequiredProperty("cromwell_jar_path");
        cromwellHost = getRequiredProperty("cromwell_host");
        wdlWorkflow = getRequiredProperty("wdl_workflow");
        wdlInputs = getRequiredProperty("wdl_inputs");
        wdlOptions = getOptionalProperty("wdl_options", null);
        wdlDepsZip = getOptionalProperty("wdl_deps_zip", null);
        outputFilesDefinition = getOptionalProperty("output_files", null);
        manualOutput = Boolean.valueOf(getRequiredProperty("manual_output"));
    }

    @Override
    public void setupDirectory() {
        //since setupDirectory is the first method run, we use it to initialize variables too.
        init();

        this.addDirectory(tmpDir);
        this.addDirectory(outputDir);
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {

        // job write wdl workflow to file
        String wdlWorkflowFile;
        Job writeWdlToFile = null;
        if (wdlWorkflow.startsWith("file://")) {
            wdlWorkflowFile = wdlWorkflow;
        } else {
            wdlWorkflowFile = tmpDir + "workflow.wdl";
            writeWdlToFile = newJob("write_wdl");
            writeWdlToFile.setMaxMemory("4000");
            writeWdlToFile.getCommand().setArguments(writeStringToFile(wdlWorkflow, wdlWorkflowFile));
        }

        // job write wdl inputs to file
        String wdlWorkflowInputsFile;
        Job writeWdlInputsToFile = null;
        if (wdlInputs.startsWith("file://")) {
            //cromwell current does not support inputs urls, remove file:// prefix
            wdlWorkflowInputsFile = wdlInputs.replaceFirst("^file://", "");
        } else {
            wdlWorkflowInputsFile = tmpDir + "inputs.json";
            writeWdlInputsToFile = newJob("write_inputs");
            writeWdlInputsToFile.setMaxMemory("4000");
            writeWdlInputsToFile.getCommand().setArguments(writeStringToFile(wdlInputs, wdlWorkflowInputsFile));
        }

        // job write wdl outputs definition to file
        Job writeOutputsToFile = null;
        String writeOutputListFile = null;
        List<OutputFile> outputFiles = new ArrayList<>();
        if (outputFilesDefinition != null) {
            throw new UnsupportedOperationException("Output files definition support is currently not implemented.");
        }

        // job execute wdl
        String workflowIdPath = tmpDir + "workflow_id";
        Job runWdlWorkflow = newJob("run_wdl");
        runWdlWorkflow.setMaxMemory("6000");
        if (writeWdlToFile != null) {
            runWdlWorkflow.addParent(writeWdlToFile);
        }
        if (writeWdlInputsToFile != null) {
            runWdlWorkflow.addParent(writeWdlInputsToFile);
        }

        Command runWdlWorkflowCommand = runWdlWorkflow.getCommand();
        runWdlWorkflowCommand.addArgument(launchAndWaitScript);
        runWdlWorkflowCommand.addArgument("--java-path");
        runWdlWorkflowCommand.addArgument(javaPath);
        runWdlWorkflowCommand.addArgument("--cromwell-jar-path");
        runWdlWorkflowCommand.addArgument(cromwellJarPath);
        runWdlWorkflowCommand.addArgument("--cromwell-host");
        runWdlWorkflowCommand.addArgument(cromwellHost);
        runWdlWorkflowCommand.addArgument("--workflow");
        runWdlWorkflowCommand.addArgument(wdlWorkflowFile);
        runWdlWorkflowCommand.addArgument("--inputs");
        runWdlWorkflowCommand.addArgument(wdlWorkflowInputsFile);
        runWdlWorkflowCommand.addArgument("--cromwell-workflow-id-path");
        runWdlWorkflowCommand.addArgument(workflowIdPath);
        if (wdlOptions != null) {
            runWdlWorkflowCommand.addArgument("--options");
            runWdlWorkflowCommand.addArgument(wdlOptions);
        }
        if (wdlDepsZip != null) {
            runWdlWorkflowCommand.addArgument("--deps-zip");
            runWdlWorkflowCommand.addArgument(wdlDepsZip);
        }

        // job provision out files
        if (writeOutputListFile != null) {
            // get output key -> real path
            // symlink real path -> outputDir/filename
            // use outputDir+filename for provision out

            //get output files and symlink them to a known location
            Job symlinkOutput = newJob("symlink_output");
            symlinkOutput.setMaxMemory("4000");
            symlinkOutput.addParent(runWdlWorkflow);

            Command symlinkOutputCommand = symlinkOutput.getCommand();
            symlinkOutputCommand.addArgument(copyOutputScript);
            symlinkOutputCommand.addArgument("--cromwell-host");
            symlinkOutputCommand.addArgument(cromwellHost);
            symlinkOutputCommand.addArgument("--cromwell-workflow-id-path");
            symlinkOutputCommand.addArgument(workflowIdPath);
            symlinkOutputCommand.addArgument("--output-dir");
            symlinkOutputCommand.addArgument(outputDir);
            symlinkOutputCommand.addArgument("--output-file-definition");
            symlinkOutputCommand.addArgument(writeOutputListFile);

            for (OutputFile f : outputFiles) {
                SqwFile sqwJsonOutputFile = createOutputFile(outputDir + f.filename, f.metatype, manualOutput);
                if (f.isOptional) {
                    sqwJsonOutputFile.setSkipIfMissing(true);
                }
                symlinkOutput.addFile(sqwJsonOutputFile);
                sqwJsonOutputFile.setType(f.metatype);
            }
        } else {
            // provision out all output files produced by cromwell workflow
            Job provisionOut = newJob("provision_out");
            provisionOut.setMaxMemory("6000");
            provisionOut.addParent(runWdlWorkflow);

            String outputPath;
            if (manualOutput) {
                outputPath = this.getMetadata_output_file_prefix() + getMetadata_output_dir() + "/";
            } else {
                outputPath = this.getMetadata_output_file_prefix() + getMetadata_output_dir() + "/" + this.getName() + "_" + this.getVersion() + "/" + this.getRandom() + "/";
            }

            Command provisionOutCommand = provisionOut.getCommand();
            provisionOutCommand.addArgument(provisionOutScript);
            provisionOutCommand.addArgument("--java-path");
            provisionOutCommand.addArgument(javaPath);
            provisionOutCommand.addArgument("--seqware-jar-path");
            provisionOutCommand.addArgument(seqwareDistributionJarPath);
            provisionOutCommand.addArgument("--workflow-run-swid");
            provisionOutCommand.addArgument("${WORKFLOW_RUN_ACCESSION}"); //exported at runtime
            provisionOutCommand.addArgument("--processing-swid");
            provisionOutCommand.addArgument("${PROCESSING_ACCESSION}"); //exported at runtime
            provisionOutCommand.addArgument("--cromwell-host");
            provisionOutCommand.addArgument(cromwellHost);
            provisionOutCommand.addArgument("--cromwell-workflow-id-path");
            provisionOutCommand.addArgument(workflowIdPath);
            provisionOutCommand.addArgument("--output-dir");
            provisionOutCommand.addArgument(outputPath);
        }

    }

    private class OutputFile {

        String metatype;
        String outputKey;
        String filename;
        Boolean isOptional = false;
        //List<String> iusSwidsToLinkFileTo;
        //signature

        public OutputFile(String outputKey, String filename, String metatype) {
            this.outputKey = outputKey;
            this.filename = filename;
            this.metatype = metatype;
        }

        public String toCsv() {
            return outputKey + "," + filename;
        }
    }

    private List<String> writeStringToFile(String fileContents, String filePath) {
        return Arrays.asList(
                "cat << END_OF_FILE_CONTENTS" + ">",
                filePath + "\n" + fileContents + "\nEND_OF_FILE_CONTENTS\n");
    }

    private String getRequiredProperty(String key) {
        String value = getProperty(key);
        if (value.trim().isEmpty()) {
            logger.error("Property = [" + key + "] is empty.");
            setWorkflowInvalid();
        }
        return value;
    }
}

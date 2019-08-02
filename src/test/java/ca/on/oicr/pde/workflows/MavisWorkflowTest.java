package ca.on.oicr.pde.workflows;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import ca.on.oicr.pde.testing.workflow.DryRun;
import ca.on.oicr.pde.testing.workflow.TestDefinition;
import net.sourceforge.seqware.common.util.maptools.MapTools;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class MavisWorkflowTest {

    public MavisWorkflowTest() {

    }

    private MavisWorkflow getWorkflow() throws IOException {
        File defaultIniFile = new File(System.getProperty("bundleDirectory") + "/config/defaults.ini");
        String defaultIniFileContents = FileUtils.readFileToString(defaultIniFile);

        MavisWorkflow wf = new MavisWorkflow();
        wf.setConfigs(MapTools.iniString2Map(defaultIniFileContents));

        return wf;
    }

    @Test
    public void validateDevelopmentTestCases() throws IllegalAccessException, InstantiationException, IOException, Exception {
        TestDefinition td = TestDefinition.buildFromJson(FileUtils.readFileToString(new File("src/test/resources/developmentTests.json")));
        for (TestDefinition.Test t : td.getTests()) {
            Map<String, String> parameters = new HashMap<>(t.getParameters());
            parameters.put("cromwell_host", "???");

            DryRun d = new DryRun(System.getProperty("bundleDirectory"), parameters, MavisWorkflow.class);
            d.buildWorkflowModel();
            d.validateWorkflow();
        }
    }

    @Test
    public void testInit() throws IOException {

        Map<String, String> config = new HashMap<>();
        config.put("wdl_inputs", "{\"key\":\"value\"}");
        config.put("cromwell_host", "???");

        MavisWorkflow wf = getWorkflow();
        wf.getConfigs().putAll(config);
        wf.setupDirectory();

        Assert.assertEquals(wf.getProperty("output_prefix"), "./");
        Assert.assertEquals(wf.getProperty("manual_output"), "false");
    }

}

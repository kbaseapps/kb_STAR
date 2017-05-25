
package us.kbase.star;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: STARalignResults</p>
 * <pre>
 * Here is the definition of the output of the function.  The output
 * can be used by other SDK modules which call your code, or the output
 * visualizations in the Narrative.  'report_name' and 'report_ref' are
 * special output fields- if defined, the Narrative can automatically
 * render your Report.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "report_name",
    "report_ref",
    "assembly_output",
    "n_initial_contigs",
    "n_contigs_removed",
    "n_contigs_remaining"
})
public class STARalignResults {

    @JsonProperty("report_name")
    private String reportName;
    @JsonProperty("report_ref")
    private String reportRef;
    @JsonProperty("assembly_output")
    private String assemblyOutput;
    @JsonProperty("n_initial_contigs")
    private Long nInitialContigs;
    @JsonProperty("n_contigs_removed")
    private Long nContigsRemoved;
    @JsonProperty("n_contigs_remaining")
    private Long nContigsRemaining;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("report_name")
    public String getReportName() {
        return reportName;
    }

    @JsonProperty("report_name")
    public void setReportName(String reportName) {
        this.reportName = reportName;
    }

    public STARalignResults withReportName(String reportName) {
        this.reportName = reportName;
        return this;
    }

    @JsonProperty("report_ref")
    public String getReportRef() {
        return reportRef;
    }

    @JsonProperty("report_ref")
    public void setReportRef(String reportRef) {
        this.reportRef = reportRef;
    }

    public STARalignResults withReportRef(String reportRef) {
        this.reportRef = reportRef;
        return this;
    }

    @JsonProperty("assembly_output")
    public String getAssemblyOutput() {
        return assemblyOutput;
    }

    @JsonProperty("assembly_output")
    public void setAssemblyOutput(String assemblyOutput) {
        this.assemblyOutput = assemblyOutput;
    }

    public STARalignResults withAssemblyOutput(String assemblyOutput) {
        this.assemblyOutput = assemblyOutput;
        return this;
    }

    @JsonProperty("n_initial_contigs")
    public Long getNInitialContigs() {
        return nInitialContigs;
    }

    @JsonProperty("n_initial_contigs")
    public void setNInitialContigs(Long nInitialContigs) {
        this.nInitialContigs = nInitialContigs;
    }

    public STARalignResults withNInitialContigs(Long nInitialContigs) {
        this.nInitialContigs = nInitialContigs;
        return this;
    }

    @JsonProperty("n_contigs_removed")
    public Long getNContigsRemoved() {
        return nContigsRemoved;
    }

    @JsonProperty("n_contigs_removed")
    public void setNContigsRemoved(Long nContigsRemoved) {
        this.nContigsRemoved = nContigsRemoved;
    }

    public STARalignResults withNContigsRemoved(Long nContigsRemoved) {
        this.nContigsRemoved = nContigsRemoved;
        return this;
    }

    @JsonProperty("n_contigs_remaining")
    public Long getNContigsRemaining() {
        return nContigsRemaining;
    }

    @JsonProperty("n_contigs_remaining")
    public void setNContigsRemaining(Long nContigsRemaining) {
        this.nContigsRemaining = nContigsRemaining;
    }

    public STARalignResults withNContigsRemaining(Long nContigsRemaining) {
        this.nContigsRemaining = nContigsRemaining;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((("STARalignResults"+" [reportName=")+ reportName)+", reportRef=")+ reportRef)+", assemblyOutput=")+ assemblyOutput)+", nInitialContigs=")+ nInitialContigs)+", nContigsRemoved=")+ nContigsRemoved)+", nContigsRemaining=")+ nContigsRemaining)+", additionalProperties=")+ additionalProperties)+"]");
    }

}

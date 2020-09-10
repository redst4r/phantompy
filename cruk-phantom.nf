#!/usr/bin/env nextflow
// checking input args
if (!params.samples){
  exit 1, "--samples not set: should be a a file of samplenames, as stored in gs://cruk-01-kallist-nextflow"
}

if (!params.outputdir){
  exit 1, "--outputdir not set!"
}


params.kallisto_gene_map = '/home/mstrasse/resources/transcripts_to_genes.txt'
Channel
        .fromPath(params.kallisto_gene_map)
        .ifEmpty { exit 1, "kallisto gene map not found: ${params.kallisto_gene_map}" }
        .set {kallisto_gene_map}
// samples = params.samples.split()

// read the samplenames from the file, line by line
Channel
  .fromPath(params.samples)
  .splitText() {it.strip()}
  .filter(){ it != ""}
  .map {it.split()}   // splitting into sampleID and gs-path
  .set{samplelist1}


// since the splitting down there doesnt work split the samplelist and download twice
Channel
  .fromPath(params.samples)
  .splitText() {it.strip()}
  .filter(){ it != ""}
  .map {it.split()}   // splitting into sampleID and gs-path
  .set{samplelist2}

process download_kallisto_bus {
  input:
  val samplename_gs from samplelist1

  output:
  file "${samplename_gs[0]}_sorted.bus" into busfiles

  script:
  sampleid=  samplename_gs[0]
  gslocation=samplename_gs[1]
  """
  gsutil -m cp ${gslocation}/kallisto/sort_bus/bus_output/output.corrected.sort.bus ./${sampleid}_sorted.bus
  """
}


//////////////////////////////
// forking this channel does something funny, such that we cant use '*.bus' in  suspicious_CUGs
//////////////////////////////
//
// Channel
//  .from(busfiles)
//  .into{busfiles1; busfiles2}
// busfiles1.view{ "Foo emit: " + it }
// busfiles2.view{ "Bar emit: " + it }


// (busfiles1, busfiles2) = Channel.from( busfiles).into(2)
// busfiles1.subscribe{ "Foo emit: " + it }
// busfiles2.subscribe{ "Bar emit: " + it }


busfiles1 = busfiles
process suspicious_CUGs {
  publishDir "${params.outputdir}", mode: 'copy'

  input:
  file '*_sorted.bus' from busfiles1.collect()

  output:
  file 'suspicious.pkl' into suspicious
  file 'suspicious.log'


  script:
  """
  python /home/mstrasse/phantompy/cruk-phantom-cli.py suspicious --outfile suspicious.pkl *_sorted.bus > suspicious.log
  """

}
//
// process download_kallisto_bus_and_filter {
//   // need to redownload to get anonther version of the busfiles to filter
//   publishDir "${params.outputdir}", mode: 'copy'
//
//   input:
//   file 'suspicious.pkl' from suspicious
//   val samplename from samplelist2
//
//   output:
//   // this will emit multiple items!
//   file "${samplename}_sorted_filtered.bus" into filtered_bus
//
//   script:
//   """
//   gsutil -m cp gs://cruk-01-kallisto-nextflow/${samplename}/kallisto/sort_bus/bus_output/output.corrected.sort.bus ./${samplename}_sorted.bus
//   python /home/mstrasse/phantompy/cruk-phantom-cli.py filter --suspicious suspicious.pkl --inbus ./${samplename}_sorted.bus --outbus  ./${samplename}_sorted_filtered.bus
//   """
// }


process download_kallisto_bus_and_filter_v2 {
  // need to redownload to get anonther version of the busfiles to filter
  publishDir "${params.outputdir}", mode: 'copy'

  input:
  file 'suspicious.pkl' from suspicious
  val samplename_gs from samplelist2
  file t2g from kallisto_gene_map.collect()

  output:
  file "${samplename_gs[0]}" into filtered_bus


  // this does it all in one
  // - download the busfile, matrix.ec, transcripts.txt
  // - filter, and create a new folder with filtered.bus, matrix.ec, transcripts.txt
  // TODO: kallisto quantification
  // bustools count -o ${bus}_eqcount/tcc -g $t2g -e ${bus}/matrix.ec -t ${bus}/transcripts.txt ${bus}/output.corrected.sort.bus
  // bustools count -o ${bus}_genecount/gene -g $t2g -e ${bus}/matrix.ec -t ${bus}/transcripts.txt --genecounts ${bus}/output.corrected.sort.bus
  script:
  sampleid=  samplename_gs[0]
  gslocation=samplename_gs[1]
  """
  gsutil -m cp '${gslocation}/kallisto/sort_bus/bus_output/*' .
  EQCOUNT=${sampleid}/kallisto/bustools_counts/bus_output_eqcount
  GENECOUNT=${sampleid}/kallisto/bustools_counts/bus_output_genecount
  SORT=${sampleid}/kallisto/sort_bus/bus_output
  mkdir -p \$EQCOUNT
  mkdir -p \$GENECOUNT
  mkdir -p \$SORT
  python /home/mstrasse/phantompy/cruk-phantom-cli.py filter --suspicious suspicious.pkl --inbus ./output.corrected.sort.bus --outbus  \$SORT/output.corrected.sort.bus
  cp matrix.ec \$SORT
  cp transcripts.txt \$SORT

  bustools count -o \$EQCOUNT/tcc -g $t2g -e \$SORT/matrix.ec -t \$SORT/transcripts.txt \$SORT/output.corrected.sort.bus
  bustools count --genecounts  -o \$GENECOUNT/gene -g $t2g -e \$SORT/matrix.ec -t \$SORT/transcripts.txt \$SORT/output.corrected.sort.bus
  """

  // gsutil -m cp 'gs://cruk-01-kallisto-nextflow/${samplename}/kallisto/sort_bus/bus_output/*' .
  // mkdir -p bus_filtered/${samplename}
  // python /home/mstrasse/phantompy/cruk-phantom-cli.py filter --suspicious suspicious.pkl --inbus ./output.corrected.sort.bus --outbus  ./bus_filtered/${samplename}/output.corrected.sort.bus
  // cp matrix.ec ./bus_filtered/${samplename}
  // cp transcripts.txt ./bus_filtered/${samplename}
}

// process filter {
//
//   input:
//   file 'suspicious.pkl' from suspicious
//   file '*.bus' from busfiles2
//
//   output:
//   file '*.filtered.bus'
//
//   script:
//   """
//   python /home/mstrasse/phantompy/cruk-phantom-cli.py filter --suspicious suspicious.pkl --busfile $busfiles2
//   """
//
// }

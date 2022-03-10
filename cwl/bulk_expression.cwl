arguments:
- position: 0
  prefix: --out-dir
  valueFrom: output
- position: 0
  prefix: --compress-featurecounts-script
  valueFrom: /pecgs-bulk-expression/bulk_expression/shrink_featurecounts.py
- position: 0
  prefix: --generate-fpkm-script
  valueFrom: /pecgs-bulk-expression/bulk_expression/gen_fpkm.py
baseCommand:
- python
- /pecgs-bulk-expression/bulk_expression/bulk_expression.py
class: CommandLineTool
cwlVersion: v1.0
id: bulk_expression
inputs:
- id: fq_1
  inputBinding:
    itemSeparator: ','
    position: '1'
    separate: false
  type: File[]
- id: fq_2
  inputBinding:
    itemSeparator: ','
    position: '2'
    separate: false
  type: File[]
- id: star_index
  inputBinding:
    position: '0'
    prefix: --star-index
  type: Directory
- id: gtf
  inputBinding:
    position: '0'
    prefix: --gtf
  type: File
- id: gene_info
  inputBinding:
    position: '0'
    prefix: --gene-info
  type: File
- default: 16
  id: cpu
  inputBinding:
    position: '0'
    prefix: --cpu
  type: int?
- default: /miniconda/envs/bulk_expression/bin:$PATH
  id: environ_PATH
  type: string?
label: bulk_expression
outputs:
- id: output_bam
  outputBinding:
    glob: star/Aligned.sortedByCoord.out.bam
  secondaryFiles:
  - .bai
  type: File
- id: readcounts_and_fpkm_tsv
  outputBinding:
    glob: output/readcount_and_fpkm.tsv.gz
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/pecgs-bulk-expression:0.0.1
- class: ResourceRequirement
  coresMin: $(inputs.cpu)
  ramMin: 80000
- class: EnvVarRequirement
  envDef:
    PATH: $(inputs.environ_PATH)

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dir,$output);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dir:s"=>\$dir,
	"output:s"=>\$output
			) or &USAGE;
&USAGE unless ($dir and $output);
mkdir $output if (!-d $output);
mkdir "$output/image" if (!-d "$output/image");

my $dpng = "$output/image";
open Out,">$output/project.md";

#################参考基因组###########################
open In,"$dir/03.reference/ref.stat";
my  $ref= "------------------------------- ---------------------------------------\n";

while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my @line=split(/\t/,$_);
	$ref.="  $line[0]                        $line[1]                               \n";
}
close In;
$ref.="------------------------------- ---------------------------------------\n";
##########################################################################
mkdir "$dpng/pipeline" if (!-d "$dpng/pipeline");
`cp $Bin/img/pipeline1.png $dpng/pipeline/pipeline1.png`;
`cp $Bin/img/pipeline2.png $dpng/pipeline/pipeline2.png`;
`cp $Bin/img/pipeline3.png $dpng/pipeline/pipeline3.png`;

######################################################


######################数据量##########################
mkdir "$dpng/QC/" if (!-d "$dpng/QC");

`cp $dir/02.cleanFQ/QC/*.png $dpng/QC/`;

my @stat=split(/\n/,`cat $dir/02.cleanFQ/QC/*.stat`);
my %stat;
my $outline.="";
foreach my $l (@stat) {
	next if ($l =~ /^#/);
	my @info=split(/\s+/,$l);
	$stat{num}++;
	$stat{base}+=$info[8];
	$stat{Q30}+=$info[10]*$info[8];
	$stat{GC}+=$info[11]*$info[8];
	$info[11]=int($info[11]*10000)/100;
	$info[10]=int($info[10]*10000)/100;
	$outline.="|".join("|",$info[0],$info[7],$info[8],$info[11],$info[10])."|\n"
}
$stat{Q30}=int($stat{Q30}/$stat{base}*10000)/100;
$stat{GC}=int($stat{GC}/$stat{base}*10000)/100;
$stat{base}=int($stat{base}/1024/1024/1024*100)/100;
$stat{Aver}=int($stat{base}/$stat{num}*100)/100;

my @basedraw=glob("$dpng/QC/*clean.base.png");
my $id=(split(/\./,basename($basedraw[0])))[0];

#########################################################################################
mkdir "$dpng/mapping" if(!-d "$dpng/mapping");
`cp $dir/08.mapstat/*.png $dpng/mapping/`;

my @resultstat=split(/\n/,`paste $dir/08.mapstat/*.result.stat`);
my @head;
my %mapstat;
my $avermapping=0;
foreach my $l (@resultstat){
  my @info=split(/\t/,$l);
  if($l =~ /type/){
    for (my $i=1;$i<@info;$i+=2){
      push @head,$info[$i];
    }
  }else{

    my $n=0;
    for (my $i=1;$i<@info;$i+=2){
      $mapstat{$head[$n]}{$info[0]}=$info[$i];
      $n++;
    }
  }
}
my $mappingstat;
foreach my $head(@head){
  $mappingstat.="|".join("|",$head,$mapstat{$head}{"mapped ratio(%)"},$mapstat{$head}{"proper ratio(%)"},$mapstat{$head}{"duplicate ratio(%)"})."|\n";
  $avermapping+=$mapstat{$head}{"mapped ratio(%)"};
}
$avermapping=int($avermapping/(scalar @head)*10000)/100;
my $coverage;
foreach my $head(@head){
  $coverage.="|".join("|",$head,$mapstat{$head}{"cover base"},$mapstat{$head}{"genome coverage(1X)"},$mapstat{$head}{"genome coverage(5X)"},$mapstat{$head}{"average depth"})."|\n";
}
##############################################################################################
open In,"$dir/18.filter/snp.sample.xls";
my $totalsnp;
my $samplesnp;
while(<In>){
  chomp;
  next if ($_ eq ""||/^$/);
  if (/total/){
    $totalsnp=(split(/\s+/,$_))[-1];
  }else{
    s/^#//g;
    my @line=split;
    $samplesnp.="|".join("|",@line)."|\n";
    if(/sampleID/){
      $samplesnp.="|--|--|--|--|--|\n";
    }
  }
}
close In;
open In,"$dir/18.filter/snp.effects.xls";
my $snpeffects;
while (<In>){
  chomp;
  next if ($_ eq ""||/^$/ || /^#/);
  my @line=split;
  $snpeffects.="|".join("|",@line)."|\n";
  if (/Type/){
      $snpeffects.="|--|--|--|\n";
  }
}
close In;
open In,"$dir/18.filter/snp.region.xls";
my $snpregion;
while (<In>){
  chomp;
  next if ($_ eq ""||/^$/ || /^#/);
  my @line=split;
  $snpregion.="|".join("|",@line)."|\n";
  if (/Type/){
    $snpregion.="|--|--|--|\n";
  }
}
close In;
##############################################################################################
open In,"$dir/18.filter/indel.sample.xls";
my $totalindel;
my $sampleindel;
while(<In>){
  chomp;
  next if ($_ eq ""||/^$/);
  if (/total/){
    $totalindel=(split(/\s+/,$_))[-1];
  }else{
    s/^#//g;
    my @line=split;
    $sampleindel.="|".join("|",@line)."|\n";
    if(/sampleID/){
      $sampleindel.="|--|--|--|--|\n";
    }
  }
}
close In;
open In,"$dir/18.filter/indel.effects.xls";
my $indeleffects;
while (<In>){
  chomp;
  next if ($_ eq ""||/^$/ || /^#/);
  my @line=split;
  $indeleffects.="|".join("|",@line)."|\n";
  if (/Type/){
      $indeleffects.="|--|--|--|\n";
  }
}
close In;
open In,"$dir/18.filter/indel.region.xls";
my $indelregion;
while (<In>){
  chomp;
  next if ($_ eq ""||/^$/ || /^#/);
  my @line=split;
  $indelregion.="|".join("|",@line)."|\n";
  if (/Type/){
    $indelregion.="|--|--|--|\n";
  }
}
close In;
###########################################################################
open In,"$dir/15.svCalling/sv.xls";
my $totalsv;
my $samplesv;
while(<In>){
  chomp;
  next if ($_ eq ""||/^$/);
  if (/total/){
    $totalsv=(split(/\s+/,$_))[-1];
  }else{
    s/^#//g;
    my @line=split;
    $samplesv.="|".join("|",@line)."|\n";
    if(/sampleID/){
      $samplesv.="|--|--|--|--|--|\n";
    }
  }
}
close In;
#####################################################################################
my @cnvline=split(/\n/,`cat $dir/17.cnvcalling/*.xls`);
my $cnvline.="|sampleID|Duplication|Deletion|\n";
$cnvline.="|--|--|--|\n";
my $totalcnv=0;
foreach my $l (@cnvline) {
	next if ($l =~ /^#/);
	my @info=split(/\s+/,$l);
	$cnvline.="|".join("|",$info[0],$info[1],$info[2])."|\n";
	$totalcnv=$info[1]+$info[2];
}


########################################################################################











my $usage=<<"USAGE";

---
title: "基因组重测序及分析结题报告"
author: "long huang"
---

\# 项目信息

\#\# 项目研究背景

全基因组重测序（Whole Genome Sequencing，WGS）是指在已知物种基因组序列（Reference genome）信息的情况下，对物种内的不同个体进行测序，发现不同个体之间的遗传变异。如单核苷酸多态性（Single Nucleotide Polymorphism，SNP）、插入缺失（Insertion-Deletion，InDel）、结构变异（Structure Variation，SV），拷贝数变异（Copy Number Variation，CNV）等。
使用变异检测的分析结果，可以进行高级的群体遗传学进行分析，例如进化树、群体结构、遗传图谱构建等。


\#\# 参考基因组基本信息

$ref

**注：**
1. Assembel Level: 参考基因组组装水平

2. SeqNum:  组装完成的总序列数

3. Total length: 组装完成的基因组总长

4. N50: 基因组的组装指标N50，N50越长，代表组装比例越完整

5. GC%: 基因组的GC含量百分比

## 分析结果概述

本项目共对$stat{num}个样品进行测序，共获得$stat{base}G数据量，平均每个样品 $stat{Aver} G数据量，碱基质量Q30比例达到 $stat{Q30} %，GC含量为 $stat{GC} %，将数据比对到参考基因组后，统计所有样品的平均比对率为 $avermapping，经生物信息学分析，共检测到 $totalsnp 个SNP， $totalindel 个InDel， $totalsv 个SV， $totalcnv 个CNV。

# 项目流程

## 全基因重测序实验流程

样品基因组DNA检测合格后，利用超声波将DNA序列片段化形成随机片段，对片段化的DNA依次进行末端修复、3′端加A、连接测序接头后，再利用磁珠吸附富集基因组长度为400bp左右的片段，经过PCR扩增形成测序文库。建好的文库先进行文库质检，质检合格的文库用Illumina HiSeq^TM^平台进行测序，测序策略为Illumina PE150，总测序读长为300bp。建库流程见图2-1。

![](image/pipeline/pipeline1.png)


 图2-1 全基因重测序实验建库流程


## 全基因重测序生信流程


在Illumina Hiseq^TM^测序数据（Raw Data）下机之后，对下机数据进行质量控制，过滤其中低质量的数据，获得高质量的数据（Clean Data）。利用BWA软件（Li *et al.,* 2009）将Clean Data比对到参考基因组序列上，获得序列的位置归属（即BAM文件）。利用GATK的Best Practices流程（McKenna *et al.*, 2010）对BAM文件进行校正，并进行SNP和Small InDel标记的检测；利用Manta（Rausch *et al.*, 2012）和GATK（McKenna *et al.*, 2010）软件进行SV和CNV的结构变异检测。利用SNPEff软件（Cingolani P, 2012）和参考基因组的基因预测信息进行变异功能注释，得到SNP、InDel的功能注释信息。整体的分析流程如图2-2所示，使用软件及版本见表2-1。

![](image/pipeline/pipeline2.png)

 图2-2 生信分析流程图

**表2-1 生信分析软件列表**

  ----------------------- ------------------------ -----------------------
  **步骤**                **应用软件**             **版本**

  原始数据质控            Fastp                    0.19.6

  参考基因组比对          BWA                      0.7.17

  SNP、InDel变异检测      GATK                     4.0.11.0

  SV变异检测              Manta                    1.6.0

  CNV变异检测             GATK                     4.0.11.0

  变异功能注释            SNPEff                    4.3

  ----------------------- ------------------------ -----------------------

**注：**相关软件下载链接如下：

1.  fastp：<https://github.com/OpenGene/fastp>

2.  BWA：<http://bio-bwa.sourceforge.net/>

3.  GATK：<https://software.broadinstitute.org/gatk/>

4.  SAMtools：<http://www.htslib.org/>

5.  Sentieon：<https://www.sentieon.com/products/>

7.  Manta: <https://github.com/Illumina/manta>

8.  SNPEff：<http://snpeff.sourceforge.net/>


# 原始数据质控和过滤

## 原始测序数据说明

为方便测序数据的分析、发布和共享，Illumina
Hiseq^TM^平台测序得到的原始图像数据经过Base
Calling转化为序列数据，得到最原始的测序数据文件。原始数据一般存储为FASTQ格式。FASTQ格式文件可记录所测读段（Reads）的碱基及其质量分数。如图3-1所示，FASTQ格式以测序读段为单位进行存储，每条Reads在FASTQ格式文件中占四行，其中第一行和第三行由文件识别标志（Sequence
Identifiers）和读段名（ID）组成（第一行以"@"开头而第三行以"+"开头；第三行中ID可以省略，但"+"不能省略），第二行为碱基序列，第四行为对应位置碱基的测序质量分数。

![](image/pipeline/pipeline3.png){width="5.53125in" height="1.6979166666666667in"}

图3-1 读段FASTQ数据格式示例

Illumina
HiSeq^TM^测序仪一个Run有2个Flowcell，一个Flowcell中包含8个Lane，其中一个Lane包含2列，每一列又包含60个Tile，每一个Tile又会种下不同的Cluster，其产生的测序文件识别标志（Sequence
Identifiers）的详细说明如表3-1所示：

表3-1 测序文件读段识别码说明

  ------------------ ----------------------------------------------------
  **标识**           **英文描述**

  E00491             Unique instrument name

  23                 Run ID

  HVT7LCCXX          Flowcell ID

  8                  Flowcell lane

  1101               Tile number within the flowcell lane

  13352              \'x\'-coordinate of the cluster within the tile

  1520               \'y\'-coordinate of the cluster within the tile

  1                  Member of a pair, 1 or 2 (paired-end or mate-pair
                     reads only)

  N                  Y if the read fails filter (read is bad), N
                     otherwise

  0                  0 when none of the control bits are on, otherwise it
                     is an even number

  CTATAC             Index sequence
  ------------------ ----------------------------------------------------

Reads的质量分数以不同的字符来表示，在Hiseq平台中，将每个字符对应的ASCII码减去33，即为对应的测序质量值。一般地，碱基质量从0-40，即对应的ASCII码为从"!"（0+33）到"I"(40+33），碱基质量越大，可信度越高。用e表示测序错误率，用Q表示Illumina
HiSeq^TM^的碱基质量值，则有下列关系：

**Q = -10**×**lg e**

表3-2 测序错误率与测序质量值对应关系简明表

  ------------------------ ----------------------- ----------------------
  **测序错误率（e）**      **测序质量值（Q）**     **对应ASCII码**

  5%                       13                      .

  1%                       20                      5

  0.1%                     30                      ?

  0.01%                    40                      I
  ------------------------ ----------------------- ----------------------

Illumina测序属于第二代测序技术，单次运行能产生数百万级的Reads，如此海量的数据无法逐个展示每条Reads的质量情况；运用统计学的方法，对所有测序Reads的每个Cycle进行碱基分布和质量波动的统计，可以从宏观上直观地反映出样本的测序质量和文库构建质量。我们针对每一个样本的原始测序数据进行测序相关质量评估，包括A/T/G/C碱基含量分布统计和碱基错误率分布统计。

## 测序碱基含量分布统计

碱基含量分布检查一般用于检测有无AT、GC分离现象。鉴于序列的随机性和碱基互补配对的原则，理论上每个测序循环上的GC含量相等、AT含量相等，且在整个测序过程基本稳定不变，呈水平线。N为测序仪无法判断的碱基类型。本项目中$id样品的碱基含量分布图如图3-2所示，反映出该样品的文库构建质量和测序质量均可满足后续分析。


![]($dpng/QC/$id.clean.base.png){width="4.418055555555555in"
height="4.418055555555555in"}

图3-2 $id 样品的碱基组成分布图

**注：**横坐标是Reads碱基坐标，坐标表示Reads上从5\'到3\'端依次碱基的排列；纵坐标是所有Reads在该测序位置A、C、G、T、N碱基分别占的百分比，不同碱基用不同的颜色表示。序列的起始位置与测序的引物接头相连，因此A、C、G、T在起始端会有所波动，后面会趋于稳定。模糊碱基N所占比例越低，说明未知碱基数越少，测序样本受系统AT偏好影响越小。虚线左侧为Read1的统计，虚线右侧为Read2的统计结果。

## 测序碱基错误率分布统计

测序错误率会随着测序序列长度的增加而升高，这是由于测序过程中化学试剂的消耗导致的，另外，由于IlluminaHiseq^TM^测序的技术特点，测序片段前端几个Cycles和末端的错误率会偏高。本项目中$id样品的测序错误率分布如图3-3所示：

![]($dpng/QC/$id.clean.qual.png){width="4.418055555555555in"
height="4.418055555555555in"}

图3-3 $id 样品的碱基错误率分布图

**注：**横坐标是Reads碱基坐标位置，表示Reads上从5\'到3\'端依次碱基的排列；纵坐标是所有Reads在该位点处碱基的平均错误率（%）。虚线左侧为双端测序的Read1的错误率分布情况，虚线右侧为Read2的错误率分布情况。

## 原始测序数据过滤

利用Illumina的建库测序平台，构建插入片段大小为400
bp左右的测序文库。按照项目合同要求进行测序，由于Illumina的原始测序数据（Raw
Data）会存在一些质量比较低的数据，所以需要进行质量过滤，获得高质量测序数据，往往过滤低质量碱基后的Reads长度会低于测序下机Reads长度，具体标准如下：

Step 1：去除Reads中的Adapter序列；

Step 2：剪切掉5'端测序质量值低于20或识别为N的碱基；

Step 3：剪切掉3'端测序质量值低于3或识别为N的碱基；

Step 4：以4个碱基为Window，剪切掉平均质量值小于20的Window中的碱基；

Step 5：去除含N的比例达到10%的Reads；

Step 6：剪切掉超过40%的碱基质量值低于15的Reads；

Step 7：舍弃去除Adapter及质量修剪后长度小于30 bp的Reads。

对质量剪切后的Clean
Data别进Reads数、总碱基数、GC含量和Q30比例的统计，详细结果见表3-3：

[原始数据质控分析结果](./image/QC/)

表3-3 测序质量统计表

  |**SampleID**|**Clean Reads**|**Clean Base**|**GC(%)**|**Q30(%)**|
  |------|-----|---------|------|---|
  $outline

**注：**

Sample ID：样本编号；

Clean Reads：高质量的Reads数；

Clean Base：原始数据过滤后剩余的高质量测序数据总碱基数；

GC(%)：Clean Data中的GC碱基占所有碱基的百分比；

Q30(%)：Clean Data中质量值大于或等于30的碱基占所有碱基的百分比。


#基因组比对

## 基因组比对效率

在本项目中，根据合同协议，我们利用BWA软件将质控后的测序片段（Clean
Reads）比对参考基因组，比对方法为MEM。表3-4为比对结果的数据统计。

表3-4 比对结果数据统计表

 |**Sample ID**|**Mapped Ratio(%)|**Properly Mapped(%)**|**Duplication Ratio(%)**|
 |------|-----|---------|------|
 $mappingstat

**注：**Sample：样品编号；

Mapped Ratio(%)：比对到基因组的Clean Reads数占所有Clean
Reads数的百分比；

Properly
Mapped(%)：双端测序序列均定位到参考基因组上且距离符合测序片段的长度的Reads数占所有Clean
Reads的百分比；

Duplication
Ratio(%)：测序数据中冗余序列（即由于PCR产生的重复Reads）的占所有Clean
Reads的百分比，该处统计结果由Picard软件的MarkDuplicate给出。

## 插入片段分布统计

通过检测双端序列在参考基因组上的起止位置，可以得到样品DNA打断后得到的测序片段的实际大小，即插入片段大小（Insert Size），是生物信息分析时的一个重要参数。插入片段大小的分布一般符合正态分布，且只有一个单峰。$id 样品的插入片段长度分布如图3-4所示，插入片段长度分布符合正态分布，中心值在[350bp]{.ul}左右，说明测序数据文库构建无异常。

![]($dpng/mapping/$id.insert.png){width="3.1493055555555554in"
height="3.1493055555555554in"}

图3-4 $id 样品的插入片段长度分布图

**注：**横坐标为Reads对应的插入片段大小，纵坐标为相应插入片段大小所对应的Reads数。

## 深度分布统计

Reads锚定到参考基因组后，可以统计其对参考基因组的覆盖情况。参考基因组上被Reads覆盖到的碱基数占基因组总长度的百分比称为基因组覆盖度；基因组覆盖度可以反映变异检测的完整性，覆盖到参考基因组的区域越多，可以检测到的变异位点也越多。碱基上覆盖的Reads数为覆盖深度。基因组的覆盖深度会影响变异检测的准确性，在覆盖深度较高的区域（非重复序列区），变异检测的准确性也越高。另外，若基因组上碱基的覆盖深度分布较均匀，也说明测序随机性较好。本项目所有样品的碱基覆盖度和平均覆盖深度统计结果见表3-5。[XX]{.ul}样品的测序深度统计如图3-5所示，基因组覆盖度如图3-6所示。

表3-5 样品覆盖深度和覆盖度统计

 
  |**Sample ID**|**Covered Bases(bp)**|**Coverage 1X(%)**|**Coverage 5X(%)**|**Average Depth(X)**|   
  |------|-----|---------|------|---|                            
  $coverage

**注：**

Sample：样品编号；

Covered Bases
(bp)：覆盖基因组的序列长度数，一般随着测序深度上升覆盖长度上升；

Coverage 1X (%)：至少有一个测序碱基覆盖的碱基占基因组长度的百分比；

Coverage 5X (%)：至少有五个测序碱基覆盖的碱基占基因组长度的百分比；

Average Depth：样品的测序平均覆盖深度。

![]($dpng/mapping/$id.depth.png){width="3.5430555555555556in"
height="3.5430555555555556in"}

图3-5 $id 样品的深度分布图

**注：**横坐标表示测序深度，图中左侧的纵坐标轴（红色）对应红色曲线，表示对应深度的位点占全基因组的百分比，图中右侧的纵坐标（蓝色）对应蓝色曲线，表示小于或等于该深度的位点占全基因组的百分比。

![]($dpng/mapping/$id.genome.coverage.png){width="4.409722222222222in"
height="3.3069444444444445in"}

**图3-6 $id 样品的染色体覆盖深度分布图**

**注：**横坐标为基因组上的碱基位置，纵坐标为染色体上对应位置的平均覆盖深度取2的对数得到的值。

图3-6表明基因组被覆盖的较均匀，说明测序随机性较好。图上深度不均一的地方可能是由于重复序列、PCR偏好性、或着丝粒部分引起的。

[基因组比对统计结果](./image/mapping/)


# SNP检测和注释

单核苷酸多态性(Single Nucleotide Polymorphism，SNP)，主要是指在基因组水平上由单个核苷酸的变异所引起的序列多态性，是基因组上多态性最高的遗传变异之一。 SNP的变异类型分为转换和颠换两种，同种类型碱基（嘌呤与嘌呤、嘧啶与嘧啶）之间的突变称为转换（Transition）；不同类型碱基（嘌呤与嘧啶）之间的突变称为颠换（Transversion）。一般转换比颠换更容易发生，所以转换/颠换（Ti/Tv）的比值一般大于1，具体比值和所测物种有关。

### SNP检测

利用GATK的Best Practices软件处理比对结果（BAM文件），利用GATK的Haplotyper方法进行SNP检测，过滤条件按照GATK推荐的参数进行，具体参见GATK官网介绍。

与参考基因组进行比对，$stat{num} 样品中共检测到 $totalsnp 个SNP。统计结果见表3-6：

表3-6 SNP数据统计表

  $samplesnp

**注：**

Sample ID：样本编号；

Total：检测到的单核苷酸多态性位点的数量，在本表中表示材料与参考基因组之间的单核苷酸差异；

Ts/Tv：转换型SNP（Transition）和颠换型SNP（Transversion）的比值；

Het：杂合分型的SNP位点总数；

Hom：纯合分型的SNP位点总数；

### SNP位置分布

采用SnpEff软件结合本项目XX的基因组注释信息对变异位点进行注释，SnpEff会根据基因组的基因和功能区域的分布进行分析，对每个SNP所在的位置和功能进行统计，并对每个变异类型的功能进行统计。由于一个基因可能会有若干种转录方式（转录本），在本表中SNP数量可能会多于上表中SNP的数量，在此处使用的SNP为至少有一个样品与其他样品（或参考基因组）不同的SNP。表3-7展示了所有的SNP位点的的位置信息：

**表3-7 SNP位置统计表**

  $snpregion

**注：**

Downstream：位于转录终止位点的下游区域的SNP个数及所占比例；

Exon：位于外显子区域的SNP个数及所占比例；

Intergenic：位于基因间区的SNP个数及所占比例；

Intron：位于内含子区域的SNP个数及所占比例；

Splice_site_acceptor：在内含子左侧的连接点区域的SNP个数及所占比例；

Splice_site_donor：在内含子右侧的连接点区域的SNP个数及所占比例；

Splice_site_region：距离外显子或内含子2
bp的剪切位点的SNP个数及所占比例；

Transcript：位于转录区域的SNP个数及所占比例；

Upstream：位于转录起始位点的上游区域的SNP个数及所占比例；

UTR_3\_prime：位于3\`UTR区域的SNP个数及所占比例；

UTR_5\_prime：位于5\`UTR区域的SNP个数及所占比例。

### SNP功能注释

对于在全基因组的SNP对蛋白质翻译的情况进行评估，可以了解每个SNP所具有的功能，由于一个基因可能会有若干种转录方式（转录本），在本表中SNP数量可能会多于上表中SNP的数量，在此处使用的SNP为至少有一个样品与其他样品（或参考基因组）不同的SNP。所有SNP位点对基因的转录本产生的功能影响如表3-8所示。

表3-8 全基因组区域SNP功能信息统计表

   $snpeffects

**注：**

3_prime_UTR_variant:在3'UTR区域的SNP个数及比例；

5_prime_UTR_premature_start_codon_gain_variant：在5'UTR区域,增加了一个启动子的SNP个数及比例；

5_prime_UTR_variant：在5'URT区域的SNP个数及比例；

Downstream_gene_variant ：在基因下游区域的SNP个数及比例；

Intergenic_region ：在基因间区的SNP的个数及比例；

Intragenic_variant ：在基因内非功能区的SNP个数及比例；

Intron_variant ：在内含子区域的SNP位点个数及所占比例；

Missense_variant 在外显子区域的错义突变的SNP位点个数及所占比例；

Non_coding_transcript_exon_variant：非编码转录本外显子变异体；

Non_coding_transcript_variant：非编码转录本变异；

Splice_acceptor_variant：在内含子左侧的连接点区域的SNP个数及所占比例；；

Splice_donor_variant：在内含子右侧的连接点区域的SNP个数及所占比例；

Splice_region_variant：距离外显子或内含子2bp的剪切位点的SNP个数及所占比例；

Start_lost：由于SNP的突变导致启动子缺失的SNP位点个数及所占比例；

Stop_gained 由于SNP的突变导致终止子获得的SNP位点个数及所占比例；

Stop_lost：由于SNP的突变导致终止子突变的SNP位点个数及所占比例；

Stop_retained_variant
：由于SNP突变导致终止子的编码的发生了改变的SNP位点个数及所占比例；

Synonymous_variant：同义突变的SNP位点个数及所占比例；

Upstream_gene_variant 在基因上游的SNP位点个数及所占比例。

[snp详细统计结果]($dir/18.filter/snp.summary.html)

# InDel检测和注释

## InDel检测

利用GATK的Best
Practices流程处理比对结果（BAM文件），利用GATK的Haplotyper方法进行InDel检测，过滤条件按照GATK推荐的参数进行，具体参见GATK官网介绍。

对项目样品进行InDel标记开发，这里的InDel指能够明确获得序列组成的InDel标记。最终共获得$totalindel 个Indel。本次分析统计结果如表3-10所示：

表3-10 InDel数据统计表

 $sampleindel

**注：**

Sample ID：样本编号；

Total Number：检测到的插入缺失变异的位点个数；

Heterozygosity Number：杂合分型的InDel的位点个数；

Homozygosity Number：纯合分型的InDel位点个数。


## InDel位置分布

采用SnpEff软件结合本项目XX的基因组注释信息，对检测到的InDel进行注释，由于一个基因可能会有多个转录本，故本表中InDel数量可能会多于上表中InDel的数量，在此处使用的InDel为至少有一个样品与其他样品（或参考基因组）不同的InDel。表3-11为注释结果统计：

**表3-11 InDel位置信息统计表**

 $indelregion

**注：**

Downstream：位于转录终止位点的下游区域的InDel个数及所占比例；

Exon：位于外显子区域的InDel个数及所占比例；

Gene：位于基因上的InDel个数及所占比例；

Intergenic：位于基因间区的InDel个数及所占比例；

Intron：位于内含子区域的InDel个数及所占比例；

Splice_site_acceptor：在内含子左侧的连接点区域的InDel个数及所占比例；

Splice_site_donor：在内含子右侧的连接点区域的InDel个数及所占比例；

Splice_siter_egion：距离外显子或内含子2bp的剪切位点的InDel个数及所占比例；

Transcript：位于转录区域的InDel个数及所占比例；

Upstream：位于转录起始位点的上游区域的InDel个数及所占比例；

UTR_3\_prime：位于3\`UTR区域的InDel个数及所占比例；

UTR_5\_prime：位于5\`UTR区域的InDel个数及所占比例。

## InDel功能注释

对于在全基因组的InDel对蛋白质翻译的情况进行评估，可以了解每个InDel所具有的功能，由于一个基因可能会有若干种转录方式（转录本），在本表中InDel数量可能会多于上表中InDel的数量，在此处使用的InDel为至少有一个样品与其他样品（或参考基因组）不同的InDel。表3-12为全基因组区域的InDel位点对蛋白翻译影响结果统计：

表3-12 全基因组区域InDel对蛋白翻译影响结果统计表

 $indeleffects

**注：**

3_prime_UTR_variant:在3'UTR区域的InDel个数及比例；

5_prime_UTR_truncation：在5'UTR区域断裂点

5_prime_UTR_variant：在5'URT区域的InDel个数及比例；

Bidirectional_gene_fusion：双向基因融合位点及比例

Conservative_inframe_deletion：对蛋白翻译影响小的碱基缺失类型的移码突变的InDel个数；

Conservative_inframe_insertion：对蛋白翻译影响小的碱基插入类型的移码突变的InDel个数；

Disruptive_inframe_deletion：严重影响蛋白翻译的碱基缺失类型的移码突变的InDel个数；

Disruptive_inframe_insertion：严重影响蛋白翻译的碱基插入类型的移码突变的InDel个数；

Downstream_gene_variant ：在基因下游区域的InDel个数及比例；

Exon_loss_variant：导致外显子缺失的InDel个数；

Frameshift_variant：导致移码突变的InDel个数；

Intergenic_region ：在基因间区的InDel的个数及比例；

Intragenic_variant ：在基因内非功能区的InDel个数及比例；

Intron_variant ：在内含子区域的InDel位点个数及所占比例；

Non_coding_transcript_exon_variant：位于外显子上且导致无法编码蛋白的InDel个数；

Non_coding_transcript_variant：导致无法编码蛋白的InDel个数

Splice_acceptor_variant：在内含子左侧的连接点区域的InDel个数及所占比例；；

Splice_donor_variant：在内含子右侧的连接点区域的InDel个数及所占比例；

Splice_region_variant：距离外显子或内含子2bp的剪切位点的InDel个数及所占比例；

Start_lost：由于SNP的突变导致启动子缺失的InDel位点个数及所占比例；

Stop_gained：由于InDel的突变导致终止子获得的InDel位点个数及所占比例；

Stop_lost：由于InDel的突变导致终止子缺失的InDel位点个数及所占比例；

Upstream_gene_variant：在基因上游区域的InDel个数及比例。

[indel详细统计结果]($dir/18.filter/indel.summary.html)

# SV检测

染色体结构变异（SV）是染色体变异的一种，类型包括：缺失（Deletion,
DEL）、插入（Insertion, INS）、易位（Translocation, BND）和重复（Duplication,
DUP）。用软件Delly检测样本的SV（表3-14）。

表3-14 SV检测结果统计表

  $samplesv

**注：**

Sample ID：样本编号；

DEL：发生序列缺失的SV变异个数；

INS：发生序列插入的SV变异个数；

BND：染色体序列发生易位的SV变异个数；

DUP：染色体序列发生重复的SV变异个数。

# CNV检测

拷贝数变异（CNV）是由基因组发生重排而导致的，一般指长度为1
kb以上的基因组大片段的拷贝数增加或者减少，主要表现为亚显微水平的缺失和重复，会影响到基因组的稳定性、基因的功能。本项目利用GATK进行CNV分析。表3-15为CNV检测统计结果:

表3-15 CNV预测结果统计表

 $cnvline

**注：**

Sample：样本编号；

Duplication：拷贝数的增加的CNV个数；

Deletion：拷贝数的减少的CNV个数。

USAGE
print Out $usage;


close Out;
`pandoc project.md  -o project.html   --toc  --template=$Bin/template/pandoc.html  -s --self-contained --resource-path=$Bin/template/ -f markdown+raw_html`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
Usage:
  Options:
  -dir		<file>	input workflow dir
  -output	<file>	output markdown dir
  -h         Help

USAGE
        print $usage;
        exit;
}

#use warnings;



################################################
########    SET INPUT/OUTPUT FILES    ##########
################################################
$time = localtime;
$time = uc($time);
$time =~ /^[A-Z]+\s+([A-Z]+)\s+\S+\s+\S+\s+(\d\d\d\d)/;
$month = $1; $year = $2;
$km             ="kegg_mass_".$month."_".$year.".txt";
$kf             ="kegg_formulas_".$month."_".$year.".txt";
$ako            ="all_compound_info_".$month."_".$year.".txt";
$outbioc        ="all_bioc_cpd_info_".$month."_".$year.".txt";
$log            ="create_compound_db_".$month."_".$year.".log";
open(LOG, ">>", $log)||die;
print LOG localtime;
qx{wget -N https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip};    qx{unzip -u pathbank_all_metabolites.csv.zip};
qx{wget -N https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip};      qx{unzip -u hmdb_metabolites.zip};
qx{wget -O PubChem_substance_text_chebi_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22chebi%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_chebi%22}'};
qx{wget -O PubChem_substance_text_biocyc_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22biocyc%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_biocyc%22}'};
qx{wget -O PubChem_substance_text_hmdb_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22hmdb%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_hmdb%22}'};
qx{wget -O PubChem_substance_text_kegg_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22kegg%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_kegg%22}'};
qx{wget -N https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl};
$keggcpd  = `wget -q -O - http://rest.kegg.jp/list/compound`;
################################################
################################################





##################################################################################
######   INPUT PATHBANK COMPOUNDS
#############################################
print "INPUT pathbank_all_metabolites.csv\n";
open(PBANK, "pathbank_all_metabolites.csv")||die "unable to open pathbank_all_metabolites.csv: $!\n";
while(<PBANK>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);
        $line=$_;
        if($line=~/\,PW_C\d+\,ELECTRON/){next;}
        $line=~/\,PW_C\d+\,\"([^\"]+)\",/;
        $fnd=$1;
        $fnd=~s/\,/\_/g;
        $line=~s/(\,PW_C\d+\,)\"([^\"]+)\",/$1$fnd\,/;
        $line=~/\,PW_C\d+\,([^\,]*),([^\,]*),([^\,]*),([^\,]*),([^\,]*),([^\,]*),([^\,]*),.*\,([^\,]*)$/;
        %ALTS=();
        $name=$1; $name=CleanName($name);
        $hmdb=$2; if($hmdb=~/(HMDB\d+)/){       $id=$1; while(length($id)<11){$id=~s/HMDB/HMDB0/;} $ALTS{$id}=1;}
        $kegg=$3; if($kegg=~/(C\d\d\d\d\d)/){                           $ALTS{$1}=1;}
        $cheb=$4; if($cheb=~/(\d+)/){           $id="CHEBI:".$cheb;     $ALTS{$id}=1;}
        $inch=$8; if($inch=~/([A-Z\-]+)/){      $id="INCHI:".$inch;     $ALTS{$id}=1;}
        $form=$7;
        foreach my $id (keys %ALTS){
                $id=~s/\s+//g;
                if($id !~ /\w/){next;}
                foreach my $id2 (keys %ALTS){$id2=~s/\s+//g; $CMPD_ALTS{$id}{$id2}++; $CMPD_ALTS{$id2}{$id}++; }
                $CMPD_NAMES{$id}{$name}++;
                $CMPD_FORMULA{$id}{$form}++;
}       }
$akc=keys %CMPD_ALTS;
$nkc=keys %CMPD_NAMES;
$fkc=keys %CMPD_FORMULA;
$mkc=keys %CMPD_MASS;
$ckc=keys %CMPD_CHARGE;
print LOG "PATHB akc $akc nkc $nkc fkc $fkc mkc $mkc ckc $ckc\n";
#############################################
######   DONE INPUT PATHBANK COMPOUNDS
##################################################################################




##################################################################################
######   INPUT HMDB COMPOUNDS
#############################################
print "INPUT hmdb_metabolites.zip\n";
open(INHMDB, "hmdb_metabolites.xml")||die "unable to open hmdb_metabolites.xml: $!\n";
while(<INHMDB>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);

        if($_=~/\<metabolite\>/i){$inmet=1; %ALTS=(); $mass=''; $name=''; $formula=''; $charge=''; next;}
        if($_=~/\<\/metabolite\>/i){
                $inmet=0;
                foreach my $id (keys %ALTS){
                        $id=~s/\s+//g;
                        if($id !~ /\w/){next;}
                        foreach my $id2 (keys %ALTS){$id2=~s/\s+//g; $CMPD_ALTS{$id}{$id2}++; $CMPD_ALTS{$id2}{$id}++;}
                        $CMPD_NAMES{$id}{$name}++;
                        $CMPD_MASS{$id}{$mass}++;
                        $CMPD_CHARGE{$id}{$charge}++;
                        $CMPD_FORMULA{$id}{$formula}++;
                        next;
        }       }
        if($inmet==1){          if($_=~/\<accession\>(HMDB\d+)\<\//i){ $id=$1; while(length($id)<11){$id=~s/HMDB/HMDB0/;} $ALTS{$id}=1; next;}
                                 if($_=~/\<inchikey\>([\w\-]+)\<\//i){ $id="INCHI:".$1; $ALTS{$id}=1; next;}
                                   if($_=~/\<kegg.id\>(C\d{5})\<\//i){ $id=$1;          $ALTS{$id}=1; next;}
                          if($_=~/\<pubchem.compound_id\>(\d+)\<\//i){ $id="CID:".$1;   $ALTS{$id}=1; next;}
                                     if($_=~/\<chebi.id\>(\d+)\<\//i){ $id="CHEBI:".$1; $ALTS{$id}=1; next;}
                               if($_=~/\<biocyc.id\>\W*(\w[A-Z\d\,\(\)\-\_]+)\<\//i){ $id=$1; $ALTS{$id}=1; next;}
                                         if($_=~/\<name\>(.*?)\<\//i){ $name=$1;    $name=CleanName($name); next;}
                             if($_=~/\<chemical.formula\>(.*?)\<\//i){ $formula=$1; next;}
               if($_=~/\<average.molecular.weight\>(\d+\.*\d*)\<\//i){ $mass=$1;    next;}
                          if($_=~/\<kind\>physiological.charge\<\//i){ $inchrg=1;   next;}
                   if($inchrg==1 && $_=~/\<value\>([\+\-]*\d+)\<\//i){ $inchrg=0;   $charge=$1; next;}
}       }
$akc=keys %CMPD_ALTS;
$nkc=keys %CMPD_NAMES;
$fkc=keys %CMPD_FORMULA;
$mkc=keys %CMPD_MASS;
$ckc=keys %CMPD_CHARGE;
print LOG "HMDB akc $akc nkc $nkc fkc $fkc mkc $mkc ckc $ckc\n";
#############################################
######   DONE INPUT HMDB COMPOUNDS   #######
#################################################################################




#################################################################################
######   INPUT PUBCHEM COMPOUNDS
#############################################
print "INPUT PubChem_substance_text_biocyc_summary.csv\n";
open(PUBBIOC, "PubChem_substance_text_biocyc_summary.csv")||die "Unable to open PubChem_substance_text_biocyc_summary.csv\n";
while(<PUBBIOC>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);
        $_=~s/\"\"//g;
        $_ =~ /^(\d+)\,(\d*)\,\"BioCyc\"\,\"\W*(\w[A-Z\d\,\(\)\-\_]+)\"\,\"/i;
        $cid="CID:".$2;
        $bid=$3; $bid =~s/\"//g; $bid =~s/^\W+//g;
        if($cid=~/\d/ && $bid=~/\w/){$CMPD_ALTS{$cid}{$bid}++; $CMPD_ALTS{$bid}{$cid}++; }
}

print "INPUT PubChem_substance_text_chebi_summary.csv\n";
open(PUBCHEB, "PubChem_substance_text_chebi_summary.csv")||die "Unable to open PubChem_substance_text_chebi_summary.csv\n";
while(<PUBCHEB>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);
        $_=~s/\"\"//g;
        $_ =~ /^(\d+)\,(\d*)\,.*(CHEBI\:\d+)/i;
        $cid="CID:".$2;
        $bid=$3;
        if($cid=~/\d/ && $bid=~/\w/){$CMPD_ALTS{$cid}{$bid}++; $CMPD_ALTS{$bid}{$cid}++; }
}

print "INPUT PubChem_substance_text_kegg_summary.csv\n";
open(PUBKEGG, "PubChem_substance_text_kegg_summary.csv")||die "Unable to open PubChem_substance_text_kegg_summary.csv\n";
while(<PUBKEGG>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);
        $_=~s/\"\"//g;
        $_ =~ /^(\d+)\,(\d*)\,\"KEGG\"\,\"(C\d\d\d\d\d)/i;
        $cid="CID:".$2;
        $bid=$3;
        if($cid=~/\d/ && $bid=~/\w/){ $CMPD_ALTS{$cid}{$bid}++; $CMPD_ALTS{$bid}{$cid}++; }
}

print "INPUT PubChem_substance_text_hmdb_summary.csv\n";
open(PUBHMDB, "PubChem_substance_text_hmdb_summary.csv")||die "Unable to open PubChem_substance_text_hmdb_summary.csv\n";
while(<PUBHMDB>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);
        $_=~s/\"\"//g;
        $_ =~ /^(\d+)\,(\d*)\,\"Human Metabolome Database \(HMDB\)\"\,\"(HMDB\d+)/i;
        $cid="CID:".$2;
        $bid=$3; if($bid !~/HMDB/){next;} while(length($bid)<11){$bid=~s/HMDB/HMDB0/;}
        if($cid=~/\d/ && $bid=~/\w/){$CMPD_ALTS{$cid}{$bid}++; $CMPD_ALTS{$bid}{$cid}++; }
}
$akc=keys %CMPD_ALTS;
$nkc=keys %CMPD_NAMES;
$fkc=keys %CMPD_FORMULA;
$mkc=keys %CMPD_MASS;
$ckc=keys %CMPD_CHARGE;
print LOG "PUBCHEM akc $akc nkc $nkc fkc $fkc mkc $mkc ckc $ckc\n";
#############################################
######   DONE INPUT PUBCHEM COMPOUNDS
##################################################################################




##################################################################################
######   INPUT CHEBI COMPOUNDS
#############################################
print "INPUT CHEBI chebi.owl\n";
open(INCHEB, "chebi.owl")||die "unable to open chebi.owl: $!\n";
while(<INCHEB>){
        $_ = uc($_);
        if($_ !~/\w/){next;}
        $_=~s/[\r\n]+//;
        if($_=~/\<OWL.CLASS.RDF.ABOUT\=.*(CHEBI).(\d+)/){%ALTS=(); $chebi=$1.":".$2; $ALTS{$chebi}=1; $mass=''; $charge=''; $form=''; $alt=''; $incheb=1; next;}
        if($incheb==1){
                if($_=~/\<CHEBI.CHARGE.*?\>([\_\+\d]+)/){       $char=$1; next;}
                if($_=~/\<CHEBI.FORMULA.*?\>([A-Z\d\.\(\)]+)/){ $form=$1; next;}
                if($_=~/\<CHEBI.MASS.*?\>([\d\.]+)/){           $mass=$1; next;}
                if($_=~/\<CHEBI.INCHIKEY.*?\>([A-Z\-]+)/){            $alt="INCHI:".$1;         $ALTS{$alt}++; next;}
                if($_=~/\<OBOINOWL.HASDBXREF.*?\>PUBCHEM.(\d+)/){     $alt="CID:".$1;           $ALTS{$alt}++; next;}
                if($_=~/\<OBOINOWL.HASDBXREF.*?\>METACYC\W+(\w[A-Z\d\,\(\)\-\_]+)/){$alt=$1;    $ALTS{$alt}++; next;}
                if($_=~/\<OBOINOWL.HASDBXREF.*?\>HMDB.(HMDB\d+)/){  $alt=$1; while(length($alt)<11){$alt=~s/HMDB/HMDB0/;} $ALTS{$alt}++; next;}
                if($_=~/\<OBOINOWL.HASDBXREF.*?\>KEGG.(C\d\d\d\d\d)/){$alt=$1;                  $ALTS{$alt}++; next;}
        }
        if($_=~/\<\/OWL.CLASS\>/){
                $incheb=0;
                foreach my $alt (keys %ALTS){
                        $CMPD_CHARGE{$alt}{$char}++;
                        $CMPD_FORMULA{$alt}{$form}++;
                        $CMPD_MASS{$alt}{$mass}++;
                        $CMPD_ALTS{$chebi}{$alt}++;
                        $CMPD_ALTS{$alt}{$chebi}++;
                }
                next;
}       }

print "INPUT CHEBI names\n";
$start=148;
$count=0;
qx{rm -f names.tsv*};
while($count <=5){
        $count++;
        qx{rm -f names.tsv};
        $file='https://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel'.$start.'/Flat_file_tab_delimited/names.tsv.gz';
        qx{wget -O names.tsv.gz $file};
        qx{gunzip -f names.tsv.gz};
        open(INCPDNM, "names.tsv")||last;
        while(<INCPDNM>){
                if($_!~/\w/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                @stuff=split("\t", $_);
                $cpd="CHEBI:".$stuff[1];
                $name=CleanName($stuff[4]);
                @PARTS = split("_", $name);
                $mp=0; $mn='';
                foreach $p (@PARTS){ $p=~/([A-Z]+)/; if(length($1)>$mp){ $mp=length($1); $mn=$1; }}
                $CPD_NODD{$cpd}{$name}=@PARTS;
                $CPD_NLEN{$cpd}{$name}=length($mn);
        }
        $kc=keys %CPD_NODD;
        print "on count $count start $start kc $kc\n";
        $start+=12;
        qx{rm -f names.tsv.gz};
}
if(keys %CPD_NODD < 1){ 
        $file='https://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel'.$start.'/Flat_file_tab_delimited/names.tsv.gz';
        print "issue with chebi names.tsv file path $file, go to chebi and get correct path and fix in script.\n"; die; 
}
foreach my $cpd (sort(keys %CPD_NODD)){
        #FIRST CHECK FOR NAMES >= 7 CHARACTER STRINGS
        #SORT BY INCREASING NON-WORDS
        $mxo=10000; $mxl=0; $goodname='';
        foreach my $name (sort{$CPD_NODD{$cpd}{$a}<=>$CPD_NODD{$cpd}{$b}} keys %{$CPD_NODD{$cpd}}){
                $odd = $CPD_NODD{$cpd}{$name};
                $len = $CPD_NLEN{$cpd}{$name};
                if($len < 7){next;}
                #NOW LOOK FOR NAME WITH LONGEST DESCRIPTOR
                if($odd < $mxo){ $mxo=$odd; }
                if($odd <= $mxo+1 ){
                        if($len > $mxl){$mxl=$len; $goodname=$name;}
                }
        }

        #GET LONGEST NAME REGARDLESS OF ODD IF NO STRING LONGER THAN 7
        if($goodname!~/\w/){
                foreach my $name (sort{$CPD_NLEN{$cpd}{$b}<=>$CPD_NLEN{$cpd}{$a}} keys %{$CPD_NLEN{$cpd}}){
                        $goodname=$name; last;
                }
        }

        #FIX STUPID LONG NON-NAMES
        if(length($goodname) > 49){
                $goodname =~ s/\d+[A-Z]\_//g;
                $goodname =~ s/^[^A-Z]+|[^A-Z]+$//g;
                while($goodname =~ /\_.{1,2}\_/){$goodname =~ s/\_.{1,2}\_/\_/; if(length($goodname)<7){last;}}
                $goodname="MOD-".$goodname;
        }

        $CMPD_NAMES{$cpd}{$goodname}++;
}



open(RHNMS, "names.tsv")||die "unable to open names.tsv: $!\n";
while(<RHNMS>){
        if($_ !~/\w/){next;}
        $_ = uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t", $_);
        $stuff[1]="CHEBI:".$stuff[1];
        if($stuff[3]=~/HMDB|PubChem|MetaCyc|ChEBI|KEGG COM|UniProt/i){
                $name=CleanName($stuff[4]);
                $CMPD_NAMES{$stuff[1]}{$name}++;
}       }
$akc=keys %CMPD_ALTS;
$nkc=keys %CMPD_NAMES;
$fkc=keys %CMPD_FORMULA;
$mkc=keys %CMPD_MASS;
$ckc=keys %CMPD_CHARGE;
print LOG "CHEBI akc $akc nkc $nkc fkc $fkc mkc $mkc ckc $ckc\n";
#############################################
######   DONE INPUT CHEBI COMPOUNDS
##################################################################################




##################################################################################
######   INPUT BIOCYC COMPOUNDS
#############################################
print "INPUT BIOCYC compounds.dat\n";
if(-s $outbioc){
        print "LOADING EXISTING BIOCYC DATA\n";
        open(INBIOC,$outbioc)||die;
        while(<INBIOC>){
                $_ = uc($_);
                $_=~s/[\r\n]+//;
                @stuff=split("\t");
                $id=$stuff[0];
                if($stuff[1]=~/\w/){ $CMPD_FORMULA{$id}{$stuff[1]}++; }
                if($stuff[2]=~/\w/){ $CMPD_CHARGE{$id}{$stuff[2]}++; }
                if($stuff[3]=~/\w/){ @names=split(";",$stuff[3]); foreach my $name (@names){$CMPD_NAMES{$id}{$name}++;}}
                if($stuff[4]=~/\w/){ @alts=split(";",$stuff[4]); foreach my $alt (@alts){ $CMPD_ALTS{$id}{$alt}++; $CMPD_ALTS{$alt}{$id}++; }}
        }
}
else{   #GET BIOCYC COMPOUND DATA
        $x = qx{ls ./BIOCYC_NF/*/compounds.dat};
        open(OUTBIOC, ">", $outbioc)||die;
        @CPD_DAT=split("\n", $x);
        $ac = @CPD_DAT;
        if($ac < 1){print "Missing Biocyc compounds.dat files. Please download BioCyc collection.\n"; die;}
        $dc=0;
        foreach my $file (@CPD_DAT){
                $dc++;
                print "on $dc of $ac file $file\n";
                open(NBCPD,$file)||die "unable to open $file: $!\n";
                while(<NBCPD>){
                        $_ = uc($_);
                        $_=~s/[\r\n]+//;
                        if($_ =~ /UNIQUE\-ID\s+\-\s+\W*(\w[A-Z\d\,\(\)\-\_]+)/){@IDS=(); $inrxn=1; $id=$1; push(@IDS,$id); $kegg=''; $cheb=''; $hmdb=''; $pubc=''; $inch=''; $name=''; $name2='';}
                        if($inrxn==1){
                                if($_ =~ /LIGAND\-CPD\s+\"(C\d+)\"/){           $alt=$1;                push(@IDS,$alt); $CMPD_ALTS{$id}{$alt}++; $CMPD_ALTS{$alt}{$id}++;}
                                if($_ =~ /DBLINKS\s+\-\s+\(CHEBI\s+\"(\d+)/){   $alt="CHEBI:".$1;       push(@IDS,$alt); $CMPD_ALTS{$id}{$alt}++; $CMPD_ALTS{$alt}{$id}++;}
                                if($_ =~ /DBLINKS\s+\-\s+\(PUBCHEM\s+\"(\d+)/){ $alt="CID:".$1;         push(@IDS,$alt); $CMPD_ALTS{$id}{$alt}++; $CMPD_ALTS{$alt}{$id}++;}
                                if($_ =~ /DBLINKS\s+\-\s+\(HMDB\s+\"(HMDB\d+)/){$alt=$1;
                                        while(length($alt)<11){$alt=~s/HMDB/HMDB0/;}                    push(@IDS,$alt); $CMPD_ALTS{$id}{$alt}++; $CMPD_ALTS{$alt}{$id}++;}
                                if($_ =~ /INCHI.KEY.*INCHIKEY\=([A-Z\-]+)/){    $alt="INCHI:".$1;       push(@IDS,$alt); $CMPD_ALTS{$id}{$alt}++; $CMPD_ALTS{$alt}{$id}++;}
                                if($_ =~ /TYPES\s+\-\s+.*?(\w\S+)/){            $name =CleanName($1); $CMPD_NAMES{$id}{$name}++;}
                                if($_ =~ /COMMON\-NAME\s+\-\s+(.*)/){           $name2=CleanName($1); $CMPD_NAMES{$id}{$name2}++;}
                        }
                        if($_=~/^\/\/\s*$/){
                                if($name=~/\w/ && $name2=~/\w/){ $name.=";".$name2; }
                                elsif($name!~/\w/ && $name=~/\w/){ $name=$name2;}
                                else{}
                                $alts=join(";",@IDS);
                                print OUTBIOC "$id\t\t\t$name\t$alts\n";
                                $inrxn=0;
        }       }       }

        print "INPUT BIOCYC metabolic-reactions.xml\n";
        $x = qx{ls ./BIOCYC_NF/*/metabolic-reactions.xml};
        open(OUTBIOC, ">>", $outbioc)||die;
        @MRX_XML=split("\n", $x);
        $ac = @MRX_XML;
        $dc=0;
        foreach my $file (@MRX_XML){
                $dc++;
                print "on $dc of $ac file $file\n";
                open(INXML, $file)||die "unable to open $file: $!\n";
                $incmpd=0;
                while(<INXML>){
                        if($_ !~ /\w/){next;}
                        $_ = uc($_);

                        #COMPOUND INFO
                        if($_=~/SPECIES METAID..([^\"]+)\"/i){$incmpd=1; $cpd=''; $formula=''; $charge=''; $name=''; %IDS=();}
                        if($incmpd==1){
                                if($_=~/name\=\"(A\s+|AN\s+)*([^\"]+)/i){               $name=CleanName($2);}
                                if($_=~/charge\=\"([\+\-]*\d+)\"/i){                    $charge=$1;}
                                if($_=~/\"([^\"]+)\".*?hasOnlySubstanceUnits/i){        $formula=$1;}
                                if($_=~/identifiers.org.biocyc.META.\W*(\w[A-Z\d\,\(\)\-\_]+)\"/i){     $cpd=$1;}
                                if($_=~/identifiers.org.kegg.compound.(C\d+)\"/i){                      $IDS{$1}=1;}
                                if($_=~/identifiers.org.chebi.(CHEBI.\d+)\"/i){                         $IDS{$1}=1;}
                                if($_=~/identifiers.org.inchikey.([A-Z\-]+)\"/i){$id="INCHI:".$1;       $IDS{$id}=1;}
                                if($_=~/identifiers.org.pubchem.compound.(\d+)\"/i){$id="CID:".$1;      $IDS{$id}=1;}
                                if($_=~/identifiers.org.hmdb.(HMDB\d+)\"/i){$id=$1; while(length($id)<11){$id=~s/HMDB/HMDB0/;} $IDS{$id}=1;}
                        }
                        if($_=~/\<\/SPECIES\>/i){
                                $incmpd=0;
                                if($cpd =~/\w/ && keys %IDS > 0){
                                        @ALTS=();
                                        foreach my $id (keys %IDS){
                                                push(@ALTS,$id);
                                                $CMPD_NAMES{$id}{$name}++;
                                                $CMPD_FORMULA{$id}{$formula}++;
                                                $CMPD_CHARGE{$id}{$charge}++;
                                                $CMPD_ALTS{$id}{$cpd}++;
                                                $CMPD_ALTS{$cpd}{$id}++;
                                        }
                                        $alts=join(";",@ALTS);
                                        print OUTBIOC "$cpd\t$formula\t$charge\t$name\t$alts\n";
}       }       }       }       }
$akc=keys %CMPD_ALTS;
$nkc=keys %CMPD_NAMES;
$fkc=keys %CMPD_FORMULA;
$mkc=keys %CMPD_MASS;
$ckc=keys %CMPD_CHARGE;
print LOG "BIOCYC akc $akc nkc $nkc fkc $fkc mkc $mkc ckc $ckc\n";
#############################################
######   DONE INPUT CHEBI COMPOUNDS
#####################################################################################




##################################################################################
######   INPUT KEGG COMPOUNDS
#############################################
print "GET KEGG COMPOUND NAME\n";
@KEGGCD = split("\n", $keggcpd);
%KEGGS=();
foreach my $x (@KEGGCD){
        $x=uc($x);
        (my $cpd, my $name)=split("\t",$x);
        $cpd=~s/CPD\://;
        $name=CleanName($name);
        $CMPD_NAMES{$cpd}{$name}++;
        $KEGGS{$cpd}=1;
}

print "GET KEGG COMPOUND MASS\n";
if(-s $km){
        open(INKMASS,$km)||die;
        while(<INKMASS>){
                if($_ !~/\w/){next;}
                $_=~s/[\r\n]+//;
                $_=uc($_);
                (my $cpd, my $mass)=split("\t",$_);
                $CMPD_MASS{$cpd}{$mass}++;
}       }
else{   $MW_250   = `wget -q -O - http://rest.kegg.jp/find/compound/0-250/mol_weight`;
        $MW_500   = `wget -q -O - http://rest.kegg.jp/find/compound/250-500/mol_weight`;
        $MW_1000  = `wget -q -O - http://rest.kegg.jp/find/compound/500-1000/mol_weight`;
        $MW_10000 = `wget -q -O - http://rest.kegg.jp/find/compound/1000-10000/mol_weight`;
        if($MW_250 =~ /\w/){  @Y=split("\n",$MW_250);   push(@KEGGMW,@Y);}
        if($MW_500 =~ /\w/){  @Y=split("\n",$MW_500);   push(@KEGGMW,@Y);}
        if($MW_1000 =~ /\w/){ @Y=split("\n",$MW_1000);  push(@KEGGMW,@Y);}
        if($MW_10000 =~ /\w/){@Y=split("\n",$MW_10000); push(@KEGGMW,@Y);}
        foreach my $x (@KEGGMW){ $x=uc($x); $x=~/CPD\:(C\d+)\s+(\S+)/; $CMPD_MASS{$1}{$2}++; }
        open(KMASS, ">", $km)||die "unable to open $km: $!\n";
        foreach my $cmpd (sort{$CMPD_MASS{$cmpd}{$b}<=>$CMPD_MASS{$cmpd}{$a}} keys %KEGGS){ print KMASS "$cmpd\t$CMPD_MASS{$cmpd}\n"; last;}
        undef(@KEGGMW);
}

print "GET KEGG MOLECULAR FORMULAS\n";
if(-s $kf){
        open(INKFORM,$kf)||die;
        while(<INKFORM>){
                if($_ !~/\w/){next;}
                $_=~s/[\r\n]+//;
                $_=uc($_);
                (my $cpd, my $form)=split("\t",$_);
                $CMPD_FORMULA{$cpd}{$form}++;
}       }
else{   $kfc="http://rest.kegg.jp/find/compound/";
        $kfz="/formula";
        for my $i (A..Z){        print "i $i\n";    $x=$kfc.$i.$kfz;    $y= `wget -q -O - "$x"`; @Y=split("\n",$y); push(@AF,@Y);
                for my $j (A..Z){print "ij $i$j\n"; $x=$kfc.$i.$j.$kfz; $y= `wget -q -O - "$x"`; @Y=split("\n",$y); push(@AF,@Y);}
                for my $j (a..z){print "ij $i$j\n"; $x=$kfc.$i.$j.$kfz; $y= `wget -q -O - "$x"`; @Y=split("\n",$y); push(@AF,@Y);}
                for my $j (0..9){print "ij $i$j\n"; $x=$kfc.$i.$j.$kfz; $y= `wget -q -O - "$x"`; @Y=split("\n",$y); push(@AF,@Y);}
        }
        foreach my $x (@AF){$x=uc($x); $x=~/CPD\:(C\d+)\s+(\S+)/; $CMPD_FORMULA{$1}{$2}++;}
        open(KFORM, ">", $kf)||die "unable to open $kf: $!\n";
        foreach my $cpd (sort{$CMPD_FORMULA{$cpd}{$b}<=>$CMPD_FORMULA{$cpd}{$a}} keys %KEGGS){ $CMPD_FORMULA{$cpd}=~s/\.+$//; print KFORM "$cpd\t$CMPD_FORMULA{$cpd}\n"; last;}
        undef(@AF);

        #Get Missing Formulas
        foreach my $cpd (keys %KEGGS){
                if(!exists($CMPD_FORMULA{$cpd})){
                        print "getting formula cpd $cpd\n";
                        $y= `wget -q -O - https://www.kegg.jp/entry/$cpd`;
                        @Y=split("<tr>", $y);
                        $form='';
                        foreach my $x (@Y){if($x=~/formula.*\>(.+?)\<br/i){$form = $1; $form =~ s/\.+$//; $form =~s/\s+//g;}}
                        $CMPD_FORMULA{$cpd}{$form}++;
                        print KFORM "$cpd\t$form\n";
}       }       }
$akc=keys %CMPD_ALTS;
$nkc=keys %CMPD_NAMES;
$fkc=keys %CMPD_FORMULA;
$mkc=keys %CMPD_MASS;
$ckc=keys %CMPD_CHARGE;
print LOG "KEGG akc $akc nkc $nkc fkc $fkc mkc $mkc ckc $ckc\n";
#############################################
######   DONE INPUT KEGG COMPOUNDS
#####################################################################################


$on=0;
print "Spread out formula mass and charge\n";
#SPREAD OUT FORMULA MASS AND CHARGE - MOST CPDs HAVE FORMULAS
foreach my $id (keys %CMPD_FORMULA){
        if(exists($CMPD_MASS{$id})){
                foreach my $form (keys %{$CMPD_FORMULA{$id}}){
                        if($form !~/\w/){next;}
                        foreach my $mass (keys %{$CMPD_MASS{$id}}){
                                if($mass !~/\d/){next;}
                                $FORM_MASS{$form}{$mass}++;
                                if($on%10000==0){print "id $id form $form mass $mass\n";}
        }       }       }
        if(exists($CMPD_CHARGE{$id})){
                foreach my $form (keys %{$CMPD_FORMULA{$id}}){
                        if($form !~/\w/){next;}
                        foreach my $char (keys %{$CMPD_CHARGE{$id}}){
                                if($char !~/\d/){next;}
                                $FORM_CHAR{$form}{$char}++;
                                if($on%10000==0){print "id $id form $form char $char\n";}
        }       }       }
        $on++;
}



$on=0;
print "clean up formula mass and charge\n";
#CLEAN-UP FORMULA MASSES AND CHARGES
foreach my $form (keys %FORM_MASS){
        $mass='';
        foreach my $mx ( sort{$FORM_MASS{$form}{$b}<=>$FORM_MASS{$form}{$a}} keys %{$FORM_MASS{$form}} ){ $mass=$mx; last; }
        delete($FORM_MASS{$form});
        if($on%10000==0){print "form $form mass $mass\n";} $on++;
        $FORM_MASS{$form}=$mass;
}
foreach my $form (keys %FORM_CHAR){
        $char='';
        foreach my $cx ( sort{$FORM_CHAR{$form}{$b}<=>$FORM_CHAR{$form}{$a}} keys %{$FORM_CHAR{$form}} ){ $char=$cx; last; }
        delete($FORM_CHAR{$form});
        $FORM_CHAR{$form}=$char;
}
$fckc = keys %FORM_CHAR;
$fmkc = keys %FORM_MASS;
print LOG "fckc $fckc fmkc $fmkc\n";




#SELECT THE TOP ALTS FOR EACH IDENTIFIER
$on=0;
print "loop through compound alts\n";
foreach my $id (sort(keys %CMPD_ALTS)){
        $sakc = keys %{$CMPD_ALTS{$id}};

        #FIRST SHRINK ALTS BY # REPEAT DATABASE HITS
        %SEC=(); %IDS=();
        $km=0;  $hm=0;  $bm=0;  $im=0;  $pm=0;  $cm=0;
        foreach my $alt (sort{$CMPD_ALTS{$id}{$b}<=>$CMPD_ALTS{$id}{$a}} keys %{$CMPD_ALTS{$id}}){
                if($alt eq $id){next;}

                #FIND MAX CNT OF CROSS-DATABASE MENTIONS
                   if($alt =~ /^HMDB\d+$/               && $CMPD_ALTS{$id}{$alt}>$hm){$hm=$CMPD_ALTS{$id}{$alt};}
                elsif($alt =~ /^C\d{5}$/                && $CMPD_ALTS{$id}{$alt}>$km){$km=$CMPD_ALTS{$id}{$alt};}
                elsif($alt =~ /^CID\:\d+$/              && $CMPD_ALTS{$id}{$alt}>$pm){$pm=$CMPD_ALTS{$id}{$alt};}
                elsif($alt =~ /^CHEBI\:\d+$/            && $CMPD_ALTS{$id}{$alt}>$cm){$cm=$CMPD_ALTS{$id}{$alt};}
                elsif($alt =~ /^INCHI\:[A-Z\-]+$/       && $CMPD_ALTS{$id}{$alt}>$im){$im=$CMPD_ALTS{$id}{$alt};}
                elsif($alt =~ /^[A-Z\d\,\(\)\-\_]+$/    && $CMPD_ALTS{$id}{$alt}>$bm){$bm=$CMPD_ALTS{$id}{$alt};}
                else{   delete($CMPD_ALTS{$alt}); delete($CMPD_MASS{$alt}); delete($CMPD_NAMES{$alt});
                        delete($CMPD_FORMULA{$alt}); delete($CMPD_CHARGE{$alt}); delete($CMPD_ALTS{$id}{$alt}); next;}

                #IF COUNT IS THE MAX SCORE = BEST = KEEP
                   if($alt =~ /^HMDB\d+$/){             if($CMPD_ALTS{$id}{$alt}>=$hm/2){$IDS{$alt}=1;}else{next;}}
                elsif($alt =~ /^C\d{5}$/){              if($CMPD_ALTS{$id}{$alt}>=$km/2){$IDS{$alt}=1;}else{next;}}
                elsif($alt =~ /^CID\:\d+$/){            if($CMPD_ALTS{$id}{$alt}>=$pm/2){$IDS{$alt}=1;}else{next;}}
                elsif($alt =~ /^CHEBI\:\d+$/){          if($CMPD_ALTS{$id}{$alt}>=$cm/2){$IDS{$alt}=1;}else{next;}}
                elsif($alt =~ /^INCHI\:[A-Z\-]+$/){     if($CMPD_ALTS{$id}{$alt}>=$im/2){$IDS{$alt}=1;}else{next;}}
                elsif($alt =~ /^[A-Z\d\,\(\)\-\_]+$/){  if($CMPD_ALTS{$id}{$alt}>=$bm/2){$IDS{$alt}=1;}else{next;}}
                else{next;}

                $altkc = keys %IDS;
                if($on%10000==0){print "id $id sakc $sakc altkc $altkc alt $alt hm $hm km $km pm $pm cm $cm im $im bm $bm\n";}
        }

        #GET CPD FORM
        $form=''; %FORM=();
        foreach my $dat (keys %{$CMPD_FORMULA{$id}}){ if($dat=~/\w/){$form=$dat; $FORM{$dat}=$CMPD_FORMULA{$id}{$dat};}}
        foreach my $alt (keys %IDS){ foreach my $dat (keys %{$CMPD_FORMULA{$alt}}){ if($dat=~/\w/){$FORM{$dat}+=$CMPD_FORMULA{$alt}{$dat};}}}
        foreach my $fm (sort{$FORM{$b}<=>$FORM{$a}} keys %FORM){$form=$fm; last;}

        if($on%10000==0){print "on $on id $id form $form mass $FORM_MASS{$form} char $FORM_CHAR{$form} name $name\n";} $on++;

        #RESET CPD FORM NAME AND ALTS
        delete($CMPD_FORMULA{$id});
        delete($CMPD_ALTS{$id});
        $CMPD_FORMULA{$id}=$form;
        foreach my $alt (keys %IDS){if($alt=~/\w/){$CMPD_ALTS{$id}{$alt}=1;}}
}



#EXPAND IDS AND MERGE CONTAINMENTS
#LOOP TILL ALL CONTAINMENTS FOUND
$change=1;
while($change != 0){
        $loop++; $change=0;
        $onid=0; $totcpds=keys %CMPD_ALTS;
        foreach my $id (keys %CMPD_ALTS){

                $kc=0; $onid++;
                %ALT=(); $ALT{$id}=1;
                $CMPD_ALTS{$id}{$id}=1;
                $totalts=keys %{$CMPD_ALTS{$id}};
                foreach my $alt (keys %{$CMPD_ALTS{$id}}){$ALT{$alt}=1;}

                #REPEATEDLY LOOP THROUGH ALT IDS, COMPARE SECONDARY ALTS
                while($kc != keys %ALT){
                        $onid2=0;
                        $kc=keys %ALT;
                        foreach my $id2 (keys %ALT){

                                $onid2++;
                                %SEC=(); $SEC{$id2}=1;
                                $CMPD_ALTS{$id2}{$id2}=1;
                                foreach my $sec (keys %{$CMPD_ALTS{$id2}}){ $SEC{$sec}=1; }
                                $skc=keys %SEC; $akc=keys %ALT;

                                $cnt=0; $nn=0; foreach my $sec (keys %SEC){if(exists($ALT{$sec})){$cnt++;}else{$nn++;}}


                                if($nn == 0){ # $id2 contained in $id
                                        foreach my $sec (keys %SEC){$CMPD_ALTS{$id}{$sec}=1; $CMPD_ALTS{$id2}{$sec}=1; $ALT{$sec}=1; $SEC{$sec}=1;}
                                        foreach my $alt (keys %ALT){$CMPD_ALTS{$id}{$alt}=1; $CMPD_ALTS{$id2}{$alt}=1; $ALT{$alt}=1; $SEC{$alt}=1;}
                                }
                                elsif( $cnt > 2 && $CMPD_FORMULA{$id} eq $CMPD_FORMULA{$id2}){
                                        foreach my $sec (keys %SEC){$CMPD_ALTS{$id}{$sec}=1; $CMPD_ALTS{$id2}{$sec}=1; $ALT{$sec}=1; $SEC{$sec}=1;}
                                        foreach my $alt (keys %ALT){$CMPD_ALTS{$id}{$alt}=1; $CMPD_ALTS{$id2}{$alt}=1; $ALT{$alt}=1; $SEC{$alt}=1;}
                                }
                                elsif( $cnt > 3 && $CMPD_FORMULA{$id}=="" || $CMPD_FORMULA{$id2}==""){
                                        foreach my $sec (keys %SEC){$CMPD_ALTS{$id}{$sec}=1; $CMPD_ALTS{$id2}{$sec}=1; $ALT{$sec}=1; $SEC{$sec}=1;}
                                        foreach my $alt (keys %ALT){$CMPD_ALTS{$id}{$alt}=1; $CMPD_ALTS{$id2}{$alt}=1; $ALT{$alt}=1; $SEC{$alt}=1;}
                                }
                                else{next;}


                                if(keys %SEC != $skc || keys %ALT != $akc){
                                        $nakc=keys %ALT; $nskc=keys %SEC;
                                        print "id $id id2 $id2 kc $kc akc $akc nakc $nakc skc $skc nskc $nskc\n"; }

                                print "on1 $onid of totcpds $totcpds and on2 $onid2 of totalts $totalts id $id id2 $id2\n";
                }       }
                if(keys %{$CMPD_ALTS{$id}} != $totalts){$change++;}
        }
        print "on loop $loop change $change\n";
}



#NOW HAVE CLEANED ALTS - MERGED INDISTINGUISHABLE
#PRINT OUT CPD IDS
open(OUTCPD, ">", $ako)|| die "unable to open $ako: $!\n";
print OUTCPD "id\tform\tchar\tmass\thmdb\tinch\tbioc\tkegg\tpubc\tcheb\tname\n";
$kgc=0; $hc=0; $ic=0; $bc=0; $cc=0; $pc=0; $frmc=0; $masc=0; $crgc=0; $namc=0;
foreach my $id (keys %CMPD_ALTS){
        @KEGG=(); @CHEB=(); @INCH=(); @BIOC=(); @HMDB=(); @PUBC=(); %ALTS=();
        foreach my $alt (sort(keys %{$CMPD_ALTS{$id}})){$ALTS{$alt}=1;}

        #GET FORMULA
        $form=''; %FORMS=();
        foreach my $alt (keys %ALTS){ $dat=$CMPD_FORMULA{$alt}; if($dat=~/\w/){ $FORMS{$dat}++; }}
        foreach my $fo (sort{$FORMS{$b}<=>$FORMS{$a}} keys %FORMS){ if($fo=~/\w/){ $form=$fo; last; }}

        #GET CPD NAME
        $name=''; %NAMES=(); @NAMES=();
        foreach my $alt (keys %ALTS){ foreach my $dat (keys %{$CMPD_NAMES{$alt}}){ if($dat=~/[A-Z]{4}/){$NAMES{$dat}+=$CMPD_NAMES{$alt}{$dat};}}}
        foreach my $nm (sort{$NAMES{$b}<=>$NAMES{$a}} keys %NAMES){push(@NAMES,$nm); $mx=@NAMES; if($mx>2){last;}}
        $name=join(";",@NAMES);

        #GET CPD MASS
        $mass=''; %MASS=();
        foreach my $alt (keys %ALTS){ foreach my $dat (keys %{$CMPD_MASS{$alt}}){  if($dat=~/\d/){$MASS{$dat}+=$CMPD_MASS{$alt}{$dat}; }}}
        foreach my $ma (sort{$MASS{$b}<=>$MASS{$a}} keys %MASS){ if($ma=~/\d/){ $mass=$ma; last; }}
        if($mass eq '' && $FORM_MASS{$form}=~/\w/){$mass=$FORM_MASS{$form};}

        #GET CPD CHAR
        $char=''; %CHAR=();
        foreach my $alt (keys %ALTS){ foreach my $dat (keys %{$CMPD_CHARGE{$alt}}){ if($dat=~/\d/){$CHAR{$dat}+=$CMPD_CHARGE{$alt}{$dat}; }}}
        foreach my $ch (sort{$CHAR{$b}<=>$CHAR{$a}} keys %CHAR){ if($ch=~/\d/){ $char=$ch; last; }}
        if($char eq '' && $FORM_CHAR{$form}=~/\w/){$char=$FORM_CHAR{$form};}

        #GET CPD ALTS
        @KEGG=(); @CHEB=(); @INCH=(); @BIOC=(); @HMDB=(); @PUBC=();
        foreach my $alt (sort(keys %ALTS)){
                   if($alt =~ /^HMDB\d+$/){             push(@HMDB,$alt);}
                elsif($alt =~ /^CHEBI\:\d+$/){          push(@CHEB,$alt);}
                elsif($alt =~ /^CID\:\d+$/){            push(@PUBC,$alt);}
                elsif($alt =~ /^C\d{5}$/){              push(@KEGG,$alt);}
                elsif($alt =~ /^INCHI\:[A-Z\-]+$/){     push(@INCH,$alt);}
                elsif($alt =~ /^[A-Z\d\,\(\)\-\_]+$/){  push(@BIOC,$alt);}
                else{}
        }
        $kegg=join(";",@KEGG);
        $bioc=join(";",@BIOC);
        $cheb=join(";",@CHEB);
        $hmdb=join(";",@HMDB);
        $pubc=join(";",@PUBC);
        $inch=join(";",@INCH);

        print OUTCPD "$id\t$form\t$char\t$mass\t$hmdb\t$inch\t$bioc\t$kegg\t$pubc\t$cheb\t$name\n";
        if($kegg=~/\w/){$kgc++;}
        if($hmdb=~/\w/){$hc++;}
        if($inch=~/\w/){$ic++;}
        if($bioc=~/\w/){$bc++;}
        if($cheb=~/\w/){$cc++;}
        if($pubc=~/\w/){$pc++;}
        if($form=~/\w/){$frmc++;}
        if($mass=~/\w/){$masc++;}
        if($char=~/\w/){$crgc++;}
        if($name=~/\w/){$namc++;}
}

$tc=keys %CMPD_ALTS;
print LOG "tc $tc kc $kgc hc $hc ic $ic bc $bc cc $cc pc $pc frmc $frmc masc $masc crgc $crgc namc $namc\n";
#tc 797156 kc 1 hc 265346 ic 373465 bc 72699 cc 318564 pc 383252 frmc 625206 masc 618595 crgc 624135 namc 459125


sub CleanName{
        $name = $_[0];
        $name =~ s/(CANDIDATUS|CANDIDATUAS|CANDIDATE|VOUCHERED|UNDESCRIBED|UNSCREENED|UNKNOWN|UNCULTIVATED|UNCULTURED)\s*/\_/g;
        $name =~ s/(UNIDENTIFIED|UNCLASSIFIED|CONTAMINATION|SCREENED|UNASSIGNED|PUTATIVE|\-*LIKE)\s*/\_/g;

        #remove junk punctuation/standardize
        $name =~ s/\s+/_/g;
        $name =~ s/[^\w]+/_/g;
        $name =~ s/\_+/\_/g;
        $name =~ s/(^\_+|\_+$)//g;

        return($name);
}

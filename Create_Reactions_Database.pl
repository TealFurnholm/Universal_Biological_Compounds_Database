#use warnings;


################################################
########    SET INPUT/OUTPUT FILES    ##########
################################################
$time = localtime;
$time = uc($time);
$time =~ /^[A-Z]+\s+([A-Z]+)\s+\S+\s+\S+\s+(\d\d\d\d)/;
$month = $1; $year = $2;
$ako            ="all_compound_info_".$month."_".$year.".txt";
$binfo          ="bioc_rxn_info_".$month."_".$year.".txt";
$kd             ="kegg_rxn_dirs_".$month."_".$year.".txt";
$allrxs         ="all_reaction_info_".$month."_".$year.".txt";
$allecrx        ="all_ec_rxns_".$month."_".$year.".txt";
$allkeggrx      ="all_kegg_rxns_".$month."_".$year.".txt";
$allupids       ="all_upid_rxns_".$month."_".$year.".txt";
$outmono        ="all_mono_rxns_".$month."_".$year.".txt";
open(LOG, ">>", "create_rxn_db_".$month."_".$year.".log")||die;
print LOG "$time\n";
################################################
################################################



##################################################
###########   START GET INPUT FILES   ############
##################################################
#rhea-chebi
qx{wget -N https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv};
qx{wget -N https://ftp.expasy.org/databases/rhea/tsv/rhea2xrefs.tsv};
qx{wget -N https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv};
qx{wget -N https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz};  qx{gunzip -fk rhea-reactions.txt.gz};


#biocyc
$/="\n";
$x = qx{ls ./BIOCYC_NF/*/metabolic-reactions.xml};      @MRX_XML=split("\n", $x); $mx=@MRX_XML;
$x = qx{ls ./BIOCYC_NF/*/rxn-list};                     @RXN_LST=split("\n", $x); $rl=@RXN_LST;
$x = qx{ls ./BIOCYC_NF/*/protein-links.dat};            @PRT_LNK=split("\n", $x); $pl=@PRT_LNK;
$x = qx{ls ./BIOCYC_NF/*/enzrxns.dat};                  @ENZ_RXN=split("\n", $x); $er=@ENZ_RXN;
$x = qx{ls ./BIOCYC_NF/*/reactions.dat};                @RXN_DAT=split("\n", $x); $rd=@RXN_DAT;
print LOG "mx $mx rl $rl pl $pl er $er rd $rd\n";
if($rd < 1 || $mx < 1 || $rl < 1 || $pl < 1 || $er< 1){ print "Missing biocyc files!\n"; die;}

#kegg
$keggorg  = `wget -q -O - http://rest.kegg.jp/list/organism`;
$ko2rn    = `wget -q -O - http://rest.kegg.jp/link/rn/ko`;
$ec2rn    = `wget -q -O - http://rest.kegg.jp/link/rn/ec`;
$keggrxn  = `wget -q -O - http://rest.kegg.jp/list/reaction`;
################################################
###########   END GET INPUT FILES   ############
################################################





################################################
##########      LOAD ALL CMPDS        ##########
################################################
print "LOADING COMPOUND IDS\n";
open(CPDIN,$ako)||die "unable to open $ako:$!\n";
while(<CPDIN>){
        $_ = uc($_);
        $_=~s/[\r\n]+//;
        if($_!~/\w/){next;}
        %ALTS=();
        @stuff=split("\t", $_);
        if($stuff[6]=~/\w/){@ax=split(";",$stuff[6]); foreach my $alt (@ax){push(@ALTS,$alt);}}
        if($stuff[7]=~/\w/){@ax=split(";",$stuff[7]); foreach my $alt (@ax){push(@ALTS,$alt);}}
        if($stuff[9]=~/\w/){@ax=split(";",$stuff[9]); foreach my $alt (@ax){push(@ALTS,$alt);}}
        foreach my $a1 (keys %ALTS){foreach my $a2 (keys %ALTS){$CMPD_ALTS{$a1}{$a2}=1;}}
}
################################################################################################
################################################################################################




################################################
##########      DO BIOCYC RXNS        ##########
################################################
if(-s $binfo){
        print "Load existing BIOCYC INFO\n";
        open(INBIOC, $binfo)||die "unable to open $binfo: $1\n";
        while(<INBIOC>){
                if($_ !~/\w/){next;}
                $_ = uc($_);
                $_=~s/[\r\n]+//;
                @stuff=split("\t", $_);
                $rxn=$stuff[0];
                if($rxn=~/^TRANS-RXN/){next;}
                @DIRS =split(";",$stuff[1]);    foreach my $x (@DIRS){$RXN_DIR{$rxn}{$x}++;}
                @RITES=split(";",$stuff[2]);    foreach my $x (@RITES){$LEFT_CMPD{$rxn}{$x}++;}
                @LEFTS=split(";",$stuff[3]);    foreach my $x (@LEFTS){$RIGHT_CMPD{$rxn}{$x}++;}
                @ALTS =split(";",$stuff[4]);    foreach my $x (@ALTS){ if($x=~/^TRANS-RXN/){next;} $RXN_ALTS{$rxn}{$x}++; $RXN_ALTS{$x}{$rxn}++;}
                @IDS  =split(";",$stuff[5]);    foreach my $x (@IDS){$EC_RXN{$x}{$rxn}++;}
                @IDS  =split(";",$stuff[6]);    foreach my $x (@IDS){$UPID_RXN{$x}{$rxn}++;}
        }
        $ca=keys %CMPD_ALTS;
        $rd=keys %RXN_DIR;
        $rc=keys %RIGHT_CMPD;
        $lc=keys %LEFT_CMPD;
        $ra=keys %RXN_ALTS;
        $uc=keys %UPID_RXN;
        $er=keys %EC_RXN;
        print LOG "1.bioc-rxd rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";
}
else{
        $dc=0;
        $rd=@RXN_DAT;
        foreach my $file (@RXN_DAT){
                $dc++;
                print "on $dc of $rd file $file\n";
                open(RXDAT, $file)||die "unable to open $file:$!\n";
                %LEFT=(); %RIGHT=(); %UP=(); %IDS=();
                while(<RXDAT>){
                        $_ = uc($_);
                        $_=~s/[\r\n]+//;
                        if($_ =~ /UNIQUE-ID/){ $inrxn=1; $rxn=''; $dir=''; %LEFT=(); %RIGHT=(); %UP=(); %IDS=(); $ec='';}
                        if($inrxn==1){
                                if($_ =~ /UNIQUE\-ID\s+\-\s+[\-\+]*(\w\S+)/){$rxn=$1; if($rxn =~ /^TRANS-RXN/){ $inrxn=0; next;} $rxn=~s/^\W+//; $IDS{$rxn}=1;}
                                if($_ =~ /DBLINKS.*RHEA\s+\"(\d+)\"/){          $rhea=$1; $IDS{$rhea}=1;}
                                if($_ =~ /LIGAND\-RXN\s+\"(R\d+)\"/){           $kegg=$1; $IDS{$kegg}=1;}
                                if($_ =~ /EC\-NUMBER\s\-\sEC\-(\d+\.\d+\.\d+\.\d+)\s*$/){$ec=$1; $EC_RXN{$ec}{$rxn}++;}
                                if($_ =~ /UNIPROT\s+\"([^\"]+)\"/){$UP{$1}=1;}
                                if($_ =~ /RIGHT\s\-\s(\S+)/){$RIGHT{$1}=1;}
                                if($_ =~ /LEFT\s\-\s(\S+)/){$LEFT{$1}=1;}
                                if($_ =~ /REACTION\-DIRECTION.*(LEFT.*RIGHT)/){$dir="LTR";}
                                if($_ =~ /REACTION\-DIRECTION.*(RIGHT.*LEFT)/){$dir="RTL";}
                                if($_ =~ /REACTION\-DIRECTION.*(REVERSIBLE)/){$dir="BOTH";}
                        }
                        if($_=~/^\/\/\s*$/ && $inrxn==1){
                                $inrxn=0;
                                #correct and/or fill reactants/products and direction
                                if($dir eq "BOTH"){$RXN_DIR{$rxn}{BOTH}++;}
                                              else{$RXN_DIR{$rxn}{LTR}++;}
                                if($dir eq "RTL"){ #switch RTL so all LTR
                                        foreach my $left (keys %LEFT){  $RIGHT_CMPD{$rxn}{$left}++;}
                                        foreach my $right (keys %RIGHT){ $LEFT_CMPD{$rxn}{$right}++;}}
                                else{   foreach my $left (keys %LEFT){   $LEFT_CMPD{$rxn}{$left}++;}
                                        foreach my $right (keys %RIGHT){$RIGHT_CMPD{$rxn}{$right}++;}}
                                foreach my $upid (keys %UP){ $UPID_RXN{$upid}{$rxn}++;}
                                foreach my $xid (keys %IDS){ foreach my $yid (keys %IDS){$RXN_ALTS{$xid}{$yid}++;}}
        }       }       }
        $rd=keys %RXN_DIR;
        $rc=keys %RIGHT_CMPD;
        $lc=keys %LEFT_CMPD;
        $uc=keys %UPID_RXN;
        $ra=keys %RXN_ALTS;
        $ca=keys %CMPD_ALTS;
        $er=keys %EC_RXN;
        print LOG "1.bioc-rxd rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";

        $dc=0;
        $mx=@MRX_XML;
        foreach my $file (@MRX_XML){
                $dc++;
                print "on $dc of $mx file $file\n";
                #$file="./BIOCYC_NF/apas634452cyc/metabolic-reactions.xml";
                open(INXML, $file)||die "unable to open $file: $!\n";
                $incmpd=0;
                $inrxn=0;
                %RXNIDS=(); %CPDIDS=(); %LEFT=(); %RIGHT=(); %ID_CMPD=();
                while(<INXML>){
                        if($_ !~ /\w/){next;}
                        $_ = uc($_);

                        #COMPOUND INFO
                        if($_=~/SPECIES METAID..([^\"]+)\"/i){$chem=$1; $incmpd=1; %CPDIDS=();}
                        if($incmpd==1){
                                if($_=~/identifiers.org.biocyc.META.([^\"]+)\"/i){$CPDIDS{$1}=1;}
                                if($_=~/identifiers.org.kegg.compound.(C\d+)\"/i){$CPDIDS{$1}=1;}
                                if($_=~/identifiers.org.chebi.(CHEBI.\d+)\"/i){   $CPDIDS{$1}=1;}
                        }
                        if($_=~/\<\/SPECIES\>/i){ $incmpd=0;
                                if($chem=~/\w/ && keys %CPDIDS > 0){
                                        foreach my $id (keys %CPDIDS){$ID_CMPD{$chem}{$id}=1;}
                        }       }


                        #REACTION INFO
                        if($_=~/reaction metaid\=/i){
                                $inrxn=1; $rxn=''; %RXNIDS=(); %ECS=(); %LCPD=(); %RCPD=();
                        }
                        if($_=~/\<listOfReactants/i){$reacts=1; $prods=0;}
                        if($_=~/\<listOfProducts/i){ $reacts=0; $prods=1;}
                        if($_=~/\<\/listOfReactants/i){$reacts=0;}
                        if($_=~/\<\/listOfProducts/i){ $prods=0;}
                        if($inrxn==1){
                                if($_=~/identifiers.org.biocyc.META.([^\"]*RXN[^\"]*)\"/i){
                                        $rxn=$1; $rxn =~s/^\W+//;
                                        if($rxn =~ /^TRANS-RXN/){ $inrxn=0; $rxn='';  next;}
                                        $RXNIDS{$rxn}=1;
                                }
                                if($_=~/identifiers.org.rhea.(\d+)/i){                          $rhea=$1; $rhea=~s/^\W+//;      $RXNIDS{$rhea}=1;}
                                if($_=~/identifiers.org.kegg.reaction.(R\d+)/i){                $kegg=$1; $kegg=~s/^\W+//;      $RXNIDS{$kegg}=1;}
                                if($_=~/identifiers.org.ec.code.(\d+\.\d+\.\d+\.\d+)/i){        $ec=$1;   $ec  =~s/^\W+//;      $ECS{$ec}=1;}
                                if($reacts==1 && $_=~/speciesReference species="([^\"]+)\"/i){  $lid=$1;  $LCPD{$lid}=1; }
                                if($prods==1 && $_=~/speciesReference species="([^\"]+)\"/i){   $rid=$1;  $RCPD{$rid}=1; }
                        }

                        if($_=~/\<\/reaction/i && $inrxn==1 && $rxn=~/\w/){
                                $inrxn=0;
                                $lol=0; $rol=0; $lor=0; $ror=0; %LEFT=(); %RIGHT=();
                                foreach my $ec (keys %ECS){$EC_RXN{$ec}{$rxn}++;}
                                foreach my $lid (keys %LCPD){
                                        foreach my $left (keys %{$ID_CMPD{$lid}}){
                                                $LEFT{$left}=1;
                                                foreach my $al (keys %{$CMPD_ALTS{$left}}){
                                                        foreach my $arx (keys %{$RXN_ALTS{$rxn}}){
                                                                if(exists( $LEFT_CMPD{$arx}{$al})){$lol++;}
                                                                if(exists($RIGHT_CMPD{$arx}{$al})){$lor++;}
                                }       }       }       }
                                foreach my $rid (keys %RCPD){
                                        foreach my $right (keys %{$ID_CMPD{$rid}}){
                                                $RIGHT{$right}=1;
                                                foreach my $al (keys %{$CMPD_ALTS{$right}}){
                                                        foreach my $arx (keys %{$RXN_ALTS{$rxn}}){
                                                                if(exists( $LEFT_CMPD{$arx}{$al})){$rol++;}
                                                                if(exists($RIGHT_CMPD{$arx}{$al})){$ror++;}
                                }       }       }       }

                                if(keys %LEFT > 0 && keys %RIGHT > 0){
                                        foreach my $id (keys %RXNIDS){
                                                foreach my $id2 (keys %RXNIDS){ $RXN_ALTS{$id}{$id2}++; }
                                                if($lor > $lol && $rol >= $ror){
                                                        foreach my $left (keys %LEFT){  $RIGHT_CMPD{$id}{$left}++;}
                                                        foreach my $right (keys %RIGHT){ $LEFT_CMPD{$id}{$right}++;}
                                                }
                                                elsif($lor >= $lol && $rol > $ror){
                                                        foreach my $left (keys %LEFT){  $RIGHT_CMPD{$id}{$left}++;}
                                                        foreach my $right (keys %RIGHT){ $LEFT_CMPD{$id}{$right}++;}
                                                }
                                                else{   foreach my $left (keys %LEFT){   $LEFT_CMPD{$id}{$left}++;}
                                                        foreach my $right (keys %RIGHT){$RIGHT_CMPD{$id}{$right}++;}
        }       }       }       }       }       }

        open(OUTBIOC, ">", $binfo)||die "unable to open $binfo: $!\n";
        foreach my $up (keys %UPID_RXN){foreach my $rxn (keys %{$UPID_RXN{$up}}){ $RXN_UPID{$rxn}{$up}=$UPID_RXN{$upid}{$rxn}; }}
        foreach my $ec (keys %EC_RXN){  foreach my $rxn (keys %{$EC_RXN{$ec}}){   $RXN_EC{$rxn}{$ec}  =$EC_RXN{$ec}{$rxn}; }}
        foreach my $rxn (keys %RXN_ALTS){
                @ids=(); foreach my $idx (keys %{$RXN_DIR{$rxn}}){      push(@ids,$idx);} $dir=join(";",@ids);
                @ids=(); foreach my $idx (keys %{$RIGHT_CMPD{$rxn}}){   push(@ids,$idx);} $right=join(";",@ids);
                @ids=(); foreach my $idx (keys %{$LEFT_CMPD{$rxn}}){    push(@ids,$idx);} $left=join(";",@ids);
                @ids=(); foreach my $idx (keys %{$RXN_ALTS{$rxn}}){     push(@ids,$idx);} $alts=join(";",@ids);
                @ids=(); foreach my $idx (keys %{$RXN_EC{$rxn}}){       push(@ids,$idx);} $ecs=join(";",@ids);
                @ids=(); foreach my $idx (keys %{$RXN_UPID{$rxn}}){     push(@ids,$idx);} $ups=join(";",@ids);
                print OUTBIOC "$rxn\t$dir\t$right\t$left\t$alts\t$ecs\t$ups\n";
}       }

$rd=keys %RXN_DIR;
$rc=keys %RIGHT_CMPD;
$lc=keys %LEFT_CMPD;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$ca=keys %CMPD_ALTS;
$er=keys %EC_RXN;
print LOG "2.bioc-xml rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";
################################################################################################
################################################################################################
################################################################################################




################################################
##########        DO RHEA RXNS        ##########
################################################
print "INPUT RHEA rhea2xrefs.tsv\n";
open(RNXRF, "rhea2xrefs.tsv")|| die "unable to open rhea2xrefs.tsv: $!\n";
while(<RNXRF>){
        if($_ !~/\w/){next;}
        $_ = uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t", $_);
        $rxn=$stuff[0]; $rxn=~s/^\W+//;
        $alt=$stuff[3]; $alt=~s/^\W+//;
        if($stuff[4]=~/KEGG/ && $alt =~ /R\d\d\d\d\d/){ $RXN_ALTS{$rxn}{$alt}++; $RXN_ALTS{$alt}{$rxn}++; }
        if($stuff[4]=~/ECOCYC|METACYC|BIOCYC/){ $RXN_ALTS{$rxn}{$alt}++; $RXN_ALTS{$alt}{$rxn}++; }
        if($stuff[4]=~/^EC$/ && $stuff[3]=~/(\d+\.\d+\.\d+\.\d+)/){$EC_RXN{$1}{$rxn}++;}
}

print "INPUT RHEA rhea-directions.tsv\n";
open(INRHDIR, "rhea-directions.tsv")||die "unable to open rhea-directions.tsv: $!\n";
while(<INRHDIR>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);
        my ($main, $ltr, $rtl, $both) = split("\t",$_);
        $RHEADIR{$ltr}="LTR";
        $RHEADIR{$rtl}="RTL";
        $RXN_DIR{$main}{BOTH}++;        foreach my $alt (keys %{$RXN_ALTS{$main}}){$RXN_DIR{$alt}{BOTH}++;}
        $RXN_DIR{$ltr}{LTR}++;          foreach my $alt (keys %{$RXN_ALTS{$ltr}}){ $RXN_DIR{$alt}{LTR}++;}
        $RXN_DIR{$rtl}{RTL}++;          foreach my $alt (keys %{$RXN_ALTS{$rtl}}){ $RXN_DIR{$alt}{RTL}++;}
        $RXN_DIR{$both}{BOTH}++;        foreach my $alt (keys %{$RXN_ALTS{$both}}){$RXN_DIR{$alt}{BOTH}++;}
}

print "INPUT RHEA rhea-reactions.txt\n";
open(RHRXN, "rhea-reactions.txt")||die "unable to open rhea-reactions.txt: $!\n";
while(<RHRXN>){
        if($_ !~/\w/){next;}
        $_ = uc($_);
        $_ =~ s/[\r\n]+//;
        if($_ =~ /^ENTRY.*RHEA\D(\d+)/){$inrxn=1; $rxn=$1; $rxn=~s/^\W+//;}
        if($inrxn==1){
                if($_ =~ /^EQUATION/){
                        #Determine reaction direction based on other reaction synonyms
                        if($_ =~ /[^\<]\=\>/){$RXN_DIR{$rxn}{LTR}++;  $dir="LTR";}
                                         else{$RXN_DIR{$rxn}{BOTH}++; $dir="BOTH";}
                        if($RHEADIR{$rxn} eq "LTR"){ $dir="LTR";}
                        if($RHEADIR{$rxn} eq "RTL"){ $dir="RTL";}

                        (my $lefts, my $rights) = split("=",$_);
                        $lol=0; $rol=0; $lor=0; $ror=0;
                        @LEFTS = ( $lefts =~ /(CHEBI.\d+)/g );
                        @RIGHTS = ( $rights =~ /(CHEBI.\d+)/g );
                        if($dir eq "LTR"){
                                foreach my $left (@LEFTS){  $LEFT_CMPD{$rxn}{$left}++; }
                                foreach my $right (@RIGHTS){$RIGHT_CMPD{$rxn}{$right}++;}
                        }
                        elsif($dir eq "RTL"){
                                foreach my $right (@RIGHTS){$LEFT_CMPD{$rxn}{$right}++;}
                                foreach my $left (@LEFTS){  $RIGHT_CMPD{$rxn}{$left}++;}
                        }
                        else{
                                foreach my $left (@LEFTS){
                                        foreach my $al (keys %{$CMPD_ALTS{$left}}){
                                                foreach my $arx (keys %{$RXN_ALTS{$rxn}}){
                                                        $arx=~s/^\W+//;
                                                        if(exists($LEFT_CMPD{$arx}{$al})){$lol++;}
                                                        if(exists($RIGHT_CMPD{$arx}{$al})){$lor++;}
                                }       }       }
                                foreach my $right (@RIGHTS){
                                        foreach my $al (keys %{$CMPD_ALTS{$right}}){
                                                foreach my $arx (keys %{$RXN_ALTS{$rxn}}){
                                                        $arx=~s/^\W+//;
                                                        if(exists(  $LEFT_CMPD{$arx}{$al})){$rol++;}
                                                        if(exists(  $RIGHT_CMPD{$arx}{$al})){$ror++;}
                                }       }       }
                                #place compounds in correct orientation
                                if($lor > $lol && $rol >= $ror){
                                        foreach my $right (@RIGHTS){$LEFT_CMPD{$rxn}{$right}++;}
                                        foreach my $left (@LEFTS){  $RIGHT_CMPD{$rxn}{$left}++;}
                                }
                                elsif($lor >= $lol && $rol > $ror){
                                        foreach my $right (@RIGHTS){$LEFT_CMPD{$rxn}{$right}++;}
                                        foreach my $left (@LEFTS){  $RIGHT_CMPD{$rxn}{$left}++;}
                                }
                                else{   foreach my $left (@LEFTS){  $LEFT_CMPD{$rxn}{$left}++; }
                                        foreach my $right (@RIGHTS){$RIGHT_CMPD{$rxn}{$right}++;}
        }       }       }       }
        if($_ =~ /^\/\/\//){ $inrxn=0; $rxn=''; }
}
$rd=keys %RXN_DIR;
$rc=keys %RIGHT_CMPD;
$lc=keys %LEFT_CMPD;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$ca=keys %CMPD_ALTS;
$er=keys %EC_RXN;
print LOG "3.rhea rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";
################################################################################################
################################################################################################
################################################################################################






################################################
##########        DO KEGG RXNS        ##########
################################################
print "INPUT KEGG RXN ~ EC\n";
@EC2RN = split("\n", $ec2rn);
foreach my $x (@EC2RN){
        $x=uc($x);
        (my $ec, my $rxn)=split("\t",$x);
        $ec =~ s/EC\://;
        if($ec !~ /\d+\.\d+\.\d+\.\d+/){next;}
        $rxn =~ s/RN\://;
        $EC_RXN{$ec}{$rxn}++;
}


if(-s $kd){
        print "USING EXISTING KEGG RXN INFO $kd\n";
        open(OKRD, $kd)||die "unable to open $kd:$!\n";
        while(<OKRD>){
                #$rxn\t$lefts\t$rights\t$dirs\t$alts\n
                if($_ !~/\w/){next;}
                $_ = uc($_);
                $_ =~ s/[\r\n]+//;
                @stuff=split("\t", $_);
                $rxn=$stuff[0];
                @LEFTS=split(";",$stuff[1]);    foreach my $x (@LEFTS){$RIGHT_CMPD{$rxn}{$x}++;}
                @RITES=split(";",$stuff[2]);    foreach my $x (@RITES){$LEFT_CMPD{$rxn}{$x}++;}
                @DIRS=split(";",$stuff[3]);     foreach my $x (@DIRS){$RXN_DIR{$rxn}{$x}++;}
                @ALTS=split(";",$stuff[4]);     foreach my $x (@ALTS){$RXN_ALTS{$rxn}{$x}++; $RXN_ALTS{$x}{$rxn}++;}
}       }
else{
        print "GET KEGG REACTION DIRECTIONS\n";
        open(OKRD, ">", $kd)||die "unable to open $kd: $!\n";
        $on=0;
        @KEGGRX = split("\n", $keggrxn);
        $totkrx=@KEGGRX;
        foreach my $x (@KEGGRX){
                $on++;
                $x=uc($x);
                (my $rxn, my $junk)=split("\t",$x);
                $rxn =~ s/RN\://;
                $rxn =~s/^\W+//;


                #GET KEGG RXNS
                print "on $on of $totkrx - getting kegg $rxn rhea and l/r cpds\n";
                $re = `wget -q -O - https://www.kegg.jp/entry/$rxn`;
                $re = uc($re);
                @stuff = split("<TR>", $re);
                %LEFTS=(); %RIGHTS=();
                foreach my $line (@stuff){
                        if($line=~/rhea.db.org\D+\=(\d+)/i){$rhea=$1; $rhea=~s/\D+//g; $RXN_ALTS{$rxn}{$rhea}++; $RXN_ALTS{$rhea}{$rxn}++;}
                        if($line=~/\>EQUATION\</i){
                                (my $lline, my $rline) = split("&LT;=&GT;", $line);
                                @lefts = ( $lline =~/\>(C\d\d\d\d\d)\</g );
                                foreach my $left (@lefts){ $LEFTS{$left}=1; foreach my $al (keys %{$CMPD_ALTS{$left}}){$LEFTS{$al}=1;}}
                                @rights = ( $rline =~/\>(C\d\d\d\d\d)\</g );
                                foreach my $right (@rights){ $RIGHTS{$right}=1; foreach my $al (keys %{$CMPD_ALTS{$right}}){$RIGHTS{$al}=1;}}
                }       }

                #USE OTHER RXN ALTS TO ORIENT LEFT/RIGHT CMPDs
                $lol=0; $ror=0; $lor=0; $rol=0;
                foreach my $alt (keys %{$RXN_ALTS{$rxn}}){
                        $alt=~s/^\W+//;
                        foreach my $left (keys %{$LEFT_CMPD{$alt}}){
                                if(exists($LEFTS{$left})){ $lol++;}
                                if(exists($RIGHTS{$left})){ $lor++;}
                        }
                        foreach my $right (keys %{$RIGHT_CMPD{$alt}}){
                                if(exists($LEFTS{$right})){ $rol++;}
                                if(exists($RIGHTS{$right})){ $ror++;}
                }       }

                #place compounds in correct orientation
                   if($lor > $lol && $rol >= $ror){ foreach my $left (keys %LEFTS){ $RIGHT_CMPD{$rxn}{$left}++; } foreach my $right (keys %RIGHTS){ $LEFT_CMPD{$rxn}{$right}++; } }
                elsif($lor >= $lol && $rol > $ror){ foreach my $left (keys %LEFTS){ $RIGHT_CMPD{$rxn}{$left}++; } foreach my $right (keys %RIGHTS){ $LEFT_CMPD{$rxn}{$right}++; } }
                                              else{ foreach my $left (keys %LEFTS){  $LEFT_CMPD{$rxn}{$left}++; } foreach my $right (keys %RIGHTS){ $RIGHT_CMPD{$rxn}{$right}++;} }

                @lefts=(); @rights=(); @dirs=(); @alts=();
                foreach my $right (sort(keys %{$RIGHT_CMPD{$rxn}})){    push(@rights,$right);}  $rights =join(";",@rights);
                foreach my $left (sort(keys %{$LEFT_CMPD{$rxn}})){      push(@lefts,$left);}    $lefts  =join(";",@lefts);
                foreach my $dir (sort{$RXN_DIR{$rxn}{$b}<=>$RXN_DIR{$rxn}{$a}} keys %{$RXN_DIR{$rxn}}){ push(@dirs,$dir);} $dirs=join(";",@dirs);
                foreach my $alt (sort(keys %{$RXN_ALTS{$rxn}})){push(@alts,$alt);} $alts=join(";",@alts);
                print OKRD "$rxn\t$lefts\t$rights\t$dirs\t$alts\n";
        }
}
$rd=keys %RXN_DIR;
$rc=keys %RIGHT_CMPD;
$lc=keys %LEFT_CMPD;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$ca=keys %CMPD_ALTS;
$er=keys %EC_RXN;
print LOG "4.kegg rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";
################################################
##########       END KEGG RXNS        ##########
################################################





################################################
##########      CLEAN UP RXNS         ##########
################################################
#USE ALL COMPOUND REACTIONS TO HELP DISTRIBUTE ALL ALT REACTIONS
print "USE COMPOUND INFO TO DISTRIBUTE/EXTEND/BRIDGE BETWEEN ALT REACTIONS\n";
$rakc = keys %RXN_ALTS;
$on=0;
foreach my $id (sort(keys %RXN_ALTS)){
        $on++;
        $skc = keys %{$RXN_ALTS{$id}};

        #GET BEST OF EACH RXN TYPE
        %RXNS=(); %CPDS=(); %ALTS=(); $ALTS{$id}=1; %DIRS=();
        $ka=0;  $ra=0;  $ba=0; $dir='';
        foreach my $id2 (sort{$RXN_ALTS{$id}{$b}<=>$RXN_ALTS{$id}{$a}} keys %{$RXN_ALTS{$id}}){
                if($id eq $id2){next;}
                   if($RXN_ALTS{$id}{$id2}>$ra && $id2=~/^\d+$/){       $ra=$RXN_ALTS{$id}{$id2};}
                elsif($RXN_ALTS{$id}{$id2}>$ka && $id2=~/^R\d{5}$/){    $ka=$RXN_ALTS{$id}{$id2};}
                elsif($RXN_ALTS{$id}{$id2}>$ba && $id2=~/RXN/){         $ba=$RXN_ALTS{$id}{$id2};}
                else{   delete($RXN_ALTS{$id}{$id2}); delete($RXN_ALTS{$id2}); delete($RXN_DIR{$id2});
                        delete($LEFT_CMPD{$id2}); delete($RIGHT_CMPD{$id2}); next;}

                   if($RXN_ALTS{$id}{$id2}>=$ra/2 && $id2=~/^\d+$/){    $ALTS{$id2}=1;}
                elsif($RXN_ALTS{$id}{$id2}>=$ka/2 && $id2=~/^R\d{5}$/){ $ALTS{$id2}=1;}
                elsif($RXN_ALTS{$id}{$id2}>=$ba/2 && $id2=~/RXN/){      $ALTS{$id2}=1;}
                else{next;}
        }

        #GET ALL CPDS and CPDALTS FROM TOP RXN ALTS
        %RITE=(); %LEFT=();
        foreach my $rx (keys %ALTS){
                foreach my $dir (keys %{$RXN_DIR{$rx}}){   $DIRS{$dir}++; }
                foreach my $cpd (keys %{$LEFT_CMPD{$rx}}){ $CPDS{$cpd}++; $LEFT{$cpd}++;}
                foreach my $cpd (keys %{$RIGHT_CMPD{$rx}}){$CPDS{$cpd}++; $RITE{$cpd}++;}
        }
        foreach my $cpd (keys %CPDS){ foreach $c2 (keys %{$CMPD_ALTS{$cpd}}){ $CPDS{$c2}++; }}
        foreach my $dx (sort{$DIRS{$b}<=>$DIRS{$a}} keys %DIRS){ $dir=$dx; }


        #HAVE TOP ALTS of $id and ALL THEIR CPD ALTS,
        #LOOP THROUGH SEC ALTS - IF DOESN'T ADD NEW COMPOUNDS, THEN KEEP
        foreach my $rx (keys %ALTS){
                foreach my $alt (keys %{$RXN_ALTS{$rx}}){
                        $nn=0;
                        foreach my $cpd (keys %{$LEFT_CMPD{$alt}}){ if(!exists($CPDS{$cpd})){$nn++;}}
                        foreach my $cpd (keys %{$RIGHT_CMPD{$alt}}){if(!exists($CPDS{$cpd})){$nn++;}}
                        if($nn==0){$ALTS{$alt}=1;}
        }       }


        #UPDATE DIR, CPDS, and RXN ALTS
        delete($RXN_DIR{$id}); $RXN_DIR{$id}=$dir;
        delete($LEFT_CMPD{$id});
        delete($RIGHT_CMPD{$id});
        foreach my $cpd (keys %LEFT){ $LEFT_CMPD{$id}{$cpd}++; }
        foreach my $cpd (keys %RITE){ $RIGHT_CMPD{$id}{$cpd}++;}
        delete($RXN_ALTS{$id});
        foreach my $x (keys %ALTS){
                foreach my $y (keys %ALTS){
                        $RXN_ALTS{$x}{$y}=1;
                        $RXN_ALTS{$y}{$x}=1;
        }       }
        $akc = keys %ALTS;
        $ekc = keys %{$RXN_ALTS{$id}};
        $ckc = keys %CPDS;
        print "on $on of rakc $rakc id $id skc $skc ekc $ekc akc $akc ra $ra ba $ba ka $ka cpds $ckc\n";
}
$rd=keys %RXN_DIR;
$rc=keys %RIGHT_CMPD;
$lc=keys %LEFT_CMPD;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$ca=keys %CMPD_ALTS;
$er=keys %EC_RXN;
print LOG "5.clean1 rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";



#MERGE IDENTICAL/CONTAINED REACTIONS
#EXPAND IDS AND MERGE CONTAINMENTS
#LOOP TILL ALL CONTAINMENTS FOUND
$change=1;
while($change != 0){
        $loop++; $change=0;
        $onid=0; $totrxns=keys %RXN_ALTS;
        foreach my $id (keys %RXN_ALTS){

                $kc=0; $onid++;
                %ALT=(); $ALT{$id}=1;
                $RXN_ALTS{$id}{$id}=1;
                $totalts=keys %{$RXN_ALTS{$id}};
                foreach my $alt (keys %{$RXN_ALTS{$id}}){$ALT{$alt}=1;}

                #REPEATEDLY LOOP THROUGH ALT IDS, COMPARE SECONDARY ALTS
                while($kc != keys %ALT){
                        $onid2=0;
                        $kc=keys %ALT;
                        foreach my $id2 (keys %ALT){

                                $onid2++;
                                %SEC=(); $SEC{$id2}=1;
                                $RXN_ALTS{$id2}{$id2}=1;
                                foreach my $sec (keys %{$RXN_ALTS{$id2}}){ $SEC{$sec}=1; }
                                $skc=keys %SEC; $akc=keys %ALT;

                                $nn=0; foreach my $sec (keys %SEC){
                                if(!exists($ALT{$sec})){$nn++;}}

                                if($nn == 0){ # $id2 contained in $id
                                        foreach my $sec (keys %SEC){$RXN_ALTS{$id}{$sec}=1; $RXN_ALTS{$id2}{$sec}=1; $ALT{$sec}=1; $SEC{$sec}=1;}
                                        foreach my $alt (keys %ALT){$RXN_ALTS{$id}{$alt}=1; $RXN_ALTS{$id2}{$alt}=1; $ALT{$alt}=1; $SEC{$alt}=1;}
                                }
                                else{next;}

                                if(keys %SEC != $skc || keys %ALT != $akc){ $nakc=keys %ALT; $nskc=keys %SEC;
                                        print "id $id id2 $id2 skc $kc akc $akc nakc $nakc skc $skc nskc $nskc\n"; }

                                print "on1 $onid of totrxns $totrxns and on2 $onid2 of totalts $totalts id $id id2 $id2\n";
                }       }
                if(keys %{$RXN_ALTS{$id}} != $totalts){$change++;}
        }
        print "on loop $loop change $change\n";
}
$rd=keys %RXN_DIR;
$rc=keys %RIGHT_CMPD;
$lc=keys %LEFT_CMPD;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$ca=keys %CMPD_ALTS;
$er=keys %EC_RXN;
print LOG "6.clean2 rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";


#ORIENT LEFT/RIGHT COMPOUNDS AND REDISTRIBUTE TO REACTIONS
print "ORIENT LEFT AND RIGHT COMPOUNDS TO REACTIONS\n";
foreach my $id (sort(keys %RXN_ALTS)){

        %RXNS=(); %LEFTS=(); %RIGHTS=(); %ALL=();
        $RXNS{$id}=1;
        foreach my $id2 (keys %{$RXN_ALTS{$id}}){$RXNS{$id2}=1;}
        foreach my $id2 (keys %RXNS){ #loop through all versions/alt ids of reaction $id
                foreach my $left (keys %{$LEFT_CMPD{$id2}}){ #loop through all versions left compounds
                        $ALL{$left}=1;
                        $LEFTS{$left}++; #count each left compound
                        foreach my $alt (keys %{$CMPD_ALTS{$left}}){
                                if($alt eq $left){next;}
                                $ALL{$alt}=1; $LEFTS{$alt}++;
                        } #count each left alt compound id
                }
                foreach my $right (keys %{$RIGHT_CMPD{$id2}}){
                        $ALL{$right}=1;
                        $RIGHTS{$right}++;
                        foreach my $alt (keys %{$CMPD_ALTS{$right}}){
                                if($alt eq $right){next;}
                                $ALL{$alt}=1; $RIGHTS{$alt}++;
        }       }       }

        %NL=(); %NR=();
        foreach my $cpd (sort(keys %ALL)){
                $ll=0; $rr=0;
                foreach my $alt (keys %{$CMPD_ALTS{$cpd}}){ $ll+=$LEFTS{$alt}; }
                foreach my $alt (keys %{$CMPD_ALTS{$cpd}}){ $rr+=$RIGHTS{$alt};}

                   if($ll>$rr){$NL{$cpd}=1;}
                elsif($ll<$rr){$NR{$cpd}=1;}
                elsif($LEFT_CMPD{$id}{$cpd}>$RIGHT_CMPD{$id}{$cpd}){$NL{$cpd}=1;}
                elsif($LEFT_CMPD{$id}{$cpd}<$RIGHT_CMPD{$id}{$cpd}){$NR{$cpd}=1;}
                else{$NL{$cpd}=1; $NR{$cpd}=1;}
        }
        delete($LEFT_CMPD{$id});
        delete($RIGHT_CMPD{$id});
        foreach my $cpd (sort(keys %NL)){ $LEFT_CMPD{$id}{$cpd}=1;  }
        foreach my $cpd (sort(keys %NR)){ $RIGHT_CMPD{$id}{$cpd}=1; }
}
undef(%ALL);
undef(%CMPD_ALTS);
$rd=keys %RXN_DIR;
$rc=keys %RIGHT_CMPD;
$lc=keys %LEFT_CMPD;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$ca=keys %CMPD_ALTS;
$er=keys %EC_RXN;
print LOG "7.clean3 rd $rd rc $rc lc $lc uc $uc ra $ra ca $ca er $er\n";
################################################################################################
################################################################################################
################################################################################################






################################################
##########      PRINT OUT RXNS        ##########
################################################
print "Output reaction summary\n";
open(OUTRXN, ">", $allrxs)||die "Unable to open $allrxs: $!\n";
print OUTRXN "ID\tdir\tleft_kegg\tleft_biocyc\tleft_rhea\tright_kegg\tright_biocyc\tright_rhea\tkegg_rxn\tbiocyc_rxn\trhea_rxn\n";
foreach my $rxn (sort(keys %RXN_ALTS)){
        @KR=(); @BR=(); @RR=();
        $RXN_ALTS{$rxn}{$rxn}=1;
        foreach my $id2 (sort(keys %{$RXN_ALTS{$rxn}})){
                   if($id2=~/^\d+$/){           push(@RR,$id2);}
                elsif($id2=~/^R\d\d\d\d\d/){    push(@KR,$id2);}
                else{                           push(@BR,$id2);}
        }
        $kr=join(";",@KR); #kegg
        $br=join(";",@BR); #biocyc
        $rr=join(";",@RR); #rhea

        #get all left cpds
        @LK=(); @LB=(); @LC=();
        foreach my $cpd (sort(keys %{$RIGHT_CMPD{$rxn}})){
                         if($cpd =~ /^CHEBI/){push(@LC,$cpd);}
                elsif($cpd =~ /^C\d\d\d\d\d/){push(@LK,$cpd);}
                                         else{push(@LB,$cpd);}
        }
        $lk=join(";",@LK);
        $lb=join(";",@LB);
        $lc=join(";",@LC);

        #get all right cpds
        @RK=(); @RB=(); @RC=();
        foreach my $cpd (sort(keys %{$RIGHT_CMPD{$rxn}})){
                         if($cpd =~ /^CHEBI/){push(@RC,$cpd);}
                elsif($cpd =~ /^C\d\d\d\d\d/){push(@RK,$cpd);}
                                         else{push(@RB,$cpd);}
        }
        $rk=join(";",@RK);
        $rb=join(";",@RB);
        $rc=join(";",@RC);

        $dir=$RXN_DIR{$rxn};
        if($dir!~/\w/){$dir="BOTH";}

        print OUTRXN "$rxn\t$dir\t$lk\t$lb\t$lc\t$rk\t$rb\t$rc\t$kr\t$br\t$rr\n";
}


open(OUTECX, ">", $allecrx)||die;
foreach my $ec (keys %EC_RXN){
        if(keys %{$EC_RXN{$ec}} < 1){next;}

        #REDUCE REDUNDANT RXN ALTS
        %ALL=(); $max=0;
        foreach my $rxn (sort{$EC_RXN{$ec}{$b}<=>$EC_RXN{$ec}{$a}} keys %{$EC_RXN{$ec}}){
                if($EC_RXN{$ec}{$rxn}>$max){$max=$EC_RXN{$ec}{$rxn};}
                if($EC_RXN{$ec}{$rxn}<$max){last;}
                $isalt=0;
                foreach my $arx (keys %ALL){ if(exists($RXN_ALTS{$rxn}{$arx}) || exists($RXN_ALTS{$arx}{$rxn})){$isalt=1;}}
                if($isalt==0){$ALL{$rxn}=1;}
                if(keys %ALL > 2){last;}
        }
        @ALL=();
        foreach my $rx (keys %ALL){push(@ALL,$rx);}
        $rxn = join(";", @ALL);
        print OUTECX "$ec\t$rxn\n";
}
#########################################################
###############   END GET REACTION INFO   ###############
#########################################################






#########################################################
#########   GET UNIPROT FUNCTION CONVERSIONS   ##########
#########################################################

#GET RHEA LINKS TO UNIPROT
print "INPUT RHEA rhea2uniprot.tsv\n";
open(RHEARX, "rhea2uniprot.tsv")||die "unable to open rhea2uniprot.tsv: $!\n";
while(<RHEARX>){
        if($_ !~/\w/){next;}
        $_ = uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t", $_);
        $UPID_RXN{$stuff[3]}{$stuff[0]}++;
}
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$er=keys %EC_RXN;
print LOG "8.rhea-up uc $uc ra $ra er $er\n";


#GET KEGG GENES AND LINK THRU KOS TO RXNS
print "GET KEGG GENE KOS\n";
@Y=split("\n",$keggorg);
$count=0; $totkog=@Y;
foreach my $x (@Y){
        $count++;
        @stuff=split('\s+',$x);
        $z = `wget -q -O - http://rest.kegg.jp/link/ko/$stuff[1]`;
        $z = uc($z);
        @Z = split("\n",$z);
        foreach my $i (@Z){
                if($i !~ /K\d\d\d\d\d/){next;}
                $i=~s/\s+KO\:/\t/;
                (my $gene, my $ko)=split("\t",$i);
                $GENE_KO{$gene}{$ko}=1;
                if($count%1000==0){print "on $count of $totkog org $stuff[1] gene $gene ko $ko\n";}
}       }

print "INPUT KEGG KOs ~ RXNs \n";
@KO2RN = split("\n", $ko2rn);
foreach my $x (@KO2RN){
        $x=uc($x);
        (my $ko, my $rn)=split("\t",$x);
        $ko =~ s/KO\://;
        $rn =~ s/RN\://;
        $KO_KRXN{$ko}{$rn}=1;
}

open(OUTKGRX, ">", $allkeggrx)||die "unable to open $allkeggrx: $!\n";
foreach my $gene (keys %GENE_KO){
        @RXNS=();
        %RXNS=();
        foreach my $ko (keys %{$GENE_KO{$gene}}){
                foreach my $rxn (keys %{$KO_KRXN{$ko}}){ $RXNS{$rxn}=1; }
        }
        foreach my $rxn (keys %RXNS){ push(@RXNS, $rxn);}
        $rxns=join(";",@RXNS);
        if($rxns !~ /\w/){next;}
        print OUTKGRX "$gene\t$rxns\n";
}
$gk=keys %GENE_KO;
$kx=keys %KO_KRXN;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$er=keys %EC_RXN;
print LOG "9.kegg-gn gk $gk kx $kx uc $uc ra $ra er $er\n";
###########################################################
###########################################################



###########################################################
###########    INPUT BIOCYC PROTEIN INFO    ##############
$dc=0;
print "INPUT BIOCYC rxn-list\n";
foreach my $inrx (@RXN_LST){
        $dc++;
        print "on $dc of $ac $inrx\n";
        open(INRXL, $inrx)||die "unable to open $inrx: $!\n";
        while(<INRXL>){
                $_ = uc($_);
                if($_ !~ /\w/){next;}
                @RXL=split("KNOWN-RXN", $_);
                foreach my $rxn (@RXL){
                        @DAT=split("HIT-STRUCT", $rxn);
                        $header=shift(@DAT);
                        $header=~/\:NAME\s+(\S+)\s/;
                        $rxn=$1;
                        foreach my $dat (@DAT){
                                @HIT=split("HITDATA", $dat);
                                foreach my $hit (@HIT){
                                        if($hit =~ /QUERY\s+\"(\S+)\"/){
                                                $upid=$1;
                                                $UPID_RXN{$upid}{$rxn}++;
}       }       }       }       }       }
$gk=keys %GENE_KO;
$kx=keys %KO_KRXN;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$er=keys %EC_RXN;
print LOG "10.bioc-rxl gk $gk kx $kx uc $uc ra $ra er $er\n";



$dc=0;
$ac=@RXN_DAT;
print "INPUT BIOCYC rxn-dat\n";
foreach my $file (@RXN_DAT){
        $dc++; %UP=(); $rxn='';
        print "on $dc of $ac file $file\n";
        open(RXDAT, $file)||die "unable to open $file:$!\n";
        while(<RXDAT>){
                $_ = uc($_);
                $_=~s/[\r\n]+//;
                if($_ =~ /UNIQUE-ID/){
                        $inrxn=1;
                        $rxn='';
                        %UP=();
                }
                if($inrxn==1){
                        if($_ =~ /UNIQUE\-ID\s+\-\s+[\-\+]*(\w\S+)/){$rxn=$1;}
                        if($_ =~ /UNIPROT\s+\"([^\"]+)\"/){$UP{$1}=1;}
                }
                if($_=~/^\/\/\s*$/){
                        $inrxn=0;
                        foreach my $upid (keys %UP){ $rxn=~s/^\W+//; $UPID_RXN{$upid}{$rxn}++;}
}       }       }
$gk=keys %GENE_KO;
$kx=keys %KO_KRXN;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$er=keys %EC_RXN;
print LOG "11.bioc-rxd gk $gk kx $kx uc $uc ra $ra er $er\n";


$dc=0;
$pl=@PRT_LNK;
print "INPUT BIOCYC prt-link\n";
foreach my $inpr (@PRT_LNK){
        $dc++;
        print "on $dc of $pl $inpr\n";
        open(INPRD, $inpr)||die;
        while(<INPRD>){
                $_ = uc($_);
                @stuff=split("\t",$_);
                if($stuff[2]!~/\w/||$stuff[0]!~/MONOMER/){next;}
                $MONO_UPID{$stuff[0]}{$stuff[2]}=1;
}       }
$mp=keys %MONO_UPID;
$gk=keys %GENE_KO;
$kx=keys %KO_KRXN;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$er=keys %EC_RXN;
print LOG "12.bioc-prl mp $mp gk $gk kx $kx uc $uc ra $ra er $er\n";


print "INPUT BIOCYC enzrxns.dat\n";
$dc=0;
$on=0;
$ac=@ENZ_RXN;
foreach my $inenz (@ENZ_RXN){
        print "on $on $dc of $ac $inenz\n";
        open(INENX, $inenz)||die;
        $inrxn=0;
        $on++;
        while(<INENX>){
                $_ = uc($_);
                $_ =~ s/[\r\n]+//;
                if($_ =~ /UNIQUE-ID/){$inrxn=1; %MONO=(); %RXNS=();}
                if($inrxn==1){
                        if($_ =~ /ENZYME\s\-\s(\S*MONOMER\S*)/){$MONO{$1}=1;}
                        if($_ =~ /REACTION\s\-\s(\S+)/){$RXNS{$1}=1;}
                }
                if($_ =~ /^\/\/\s*$/){
                        $inrxn=0;
                        foreach my $mono (keys %MONO){
                                foreach my $rxn (keys %RXNS){
                                        if($rxn!~/RXN/){next;}
                                        $MONO_RXN{$mono}{$rxn}=1;
                                        foreach my $upid (keys %{$MONO_UPID{$mono}}){$UPID_RXN{$upid}{$rxn}++;}
}       }       }       }       }

        open(OUTMONO,">",$outmono)||die;
        foreach my $mono (sort(keys %MONO_RXN)){
                @RXNS=();
                foreach my $rxn (sort(keys %{$MONO_RXN{$mono}})){if($rxn=~/RXN/){push(@RXNS,$rxn);}}
                $rxns=join(";",@RXNS);
                print OUTMONO "$mono\t$rxns\n";
        }

        print "OUTPUT UPID RXNS\n";
        open(OUTUPID, ">", $allupids)||die "Unable to open $allupids: $!\n";
        foreach my $upid (keys %UPID_RXN){
                $kegg = '';
                $rhea = '';
                $bioc = '';
                foreach my $rxn (sort{$UPID_RXN{$upid}{$b}<=>$UPID_RXN{$upid}{$a}} keys %{$UPID_RXN{$upid}}){
                        if($rxn =~ /R\d\d\d\d\d$/ && $kegg eq ''){$kegg=$rxn;}
                        elsif($rxn =~ /^\d+$/ && $rhea eq ''){$rhea=$rxn;}
                        elsif($rxn =~ /R\d\d\d\d\d$/){next;}
                        elsif($rxn =~ /^\d+$/){next;}
                        elsif($rxn =~/\w/ && $bioc eq ''){$bioc=$rxn;}
                        elsif($bioc=~/\w/ && $rhea=~/\w/ && $kegg=~/\w/){last;}
                        else{next;}
                }
                print OUTUPID "$upid\t$kegg\t$rhea\t$bioc\n";
        }
$mp=keys %MONO_UPID;
$gk=keys %GENE_KO;
$kx=keys %KO_KRXN;
$uc=keys %UPID_RXN;
$ra=keys %RXN_ALTS;
$er=keys %EC_RXN;
print LOG "13.bioc-last mp $mp gk $gk kx $kx uc $uc ra $ra er $er\n";

##################################################################################################################
##################################################################################################################
##################################################################################################################


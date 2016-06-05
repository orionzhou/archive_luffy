require(dplyr)
require(ggmap)

dirw = file.path(Sys.getenv("misc2"), "urt")

## process raw data
fi = file.path(dirw, "urt_data.RData")

attach(fi)
ud = urt_data

tnames = names(ud)
yrs = sapply(strsplit(tnames, split="[_]"), "[", 1)
exs = sapply(strsplit(tnames, split="[_]"), "[", 2)
phs = sapply(strsplit(tnames, split="[_]"), "[", 3)

tp = data.frame()
phs1 = c("Height", "Lodging", "Maturity", "Oil", "Protein", "SeedQuality", "SeedSize", "YieldBuA", "YieldRank")
phs2 = c("2YearMean", "3YearMean", "4YearMean", "RegionalSummary", "Entries")
phs3 = c("Descriptive", "DescriptiveDisease")
for (i in 1:length(ud)) {
	yr = yrs[i]; ex = exs[i]; ph = phs[i]
	#DescriptiveDisease need separate processing
	if(ph %in% phs1) {
		tp1 = ud[[i]]
		if("Mean" %in% colnames(tp1)) { tp1 = tp1[,-which(colnames(tp1)=="Mean")] }
		if(colnames(tp1)[1] != 'Strain') {
			cat(tnames[i], "not start with Strain: ", colnames(tp1)[1], "\n")
			colnames(tp1)[1] = "Strain"
		}
		tp2 = reshape(tp1, direction='long', varying=list(2:ncol(tp1)), idvar='Strain', timevar="Location", v.names='Value', times=colnames(tp1)[2:ncol(tp1)])
		tp3 = cbind.data.frame(Year=yr, Test=ex, Phenotype=ph, tp2, stringsAsFactors=F)
		tp = rbind(tp, tp3)
	} else if(ph %in% phs2) {
	} else if(ph %in% phs3) {
	} else {
		cat("unknown phenotype: ", ph, "\n")
	}
}

fo = file.path(dirw, "01.tsv")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### data check
fi = file.path(dirw, "01.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti$Year = as.numeric(substr(ti$Year, 2, 5))

fo = file.path(dirw, "10.tsv")
write.table(ti, fo, sep = "\t", row.names = F, col.names = T, quote = F)

# location
locs = unique(ti$Location)
pre = sapply(locs, jj <- function(x) substr(x, start=1, stop=nchar(x)-2))
suf = sapply(locs, jj <- function(x) substr(x, start=nchar(x)-1, stop=nchar(x)))
idxs = which(suf %in% state.abb)
ts1 = data.frame(name = locs[idxs], city = pre[idxs], state = suf[idxs], stringsAsFactors = F)

idxs = which(! suf %in% state.abb)
locsc = locs[idxs]
pre = sapply(locsc, jj <- function(x) substr(x, start=1, stop=nchar(x)-3))
suf = sapply(locsc, jj <- function(x) substr(x, start=nchar(x)-2, stop=nchar(x)))
ts2 = data.frame(name = locsc, city = pre, state = suf, stringsAsFactors = F)

to = rbind(ts1, ts2)
to = to[order(to$state, to$city),]

states=to$state; cities=to$city
states[states=='QUE'] = 'QUEBEC'
lonlat = geocode(sprintf("%s %s", cities, states))
to = cbind(id=1:nrow(to), to, longitude=lonlat[,1], latitude=lonlat[,2])
fo = file.path(dirw, "11.location.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

# Test
tt = unique(ti[,c("Year", "Test")])
colnames(tt) = c("year", "name")
tt = cbind(id=1:nrow(tt), tt)
fo = file.path(dirw, "11.test.tsv")
write.table(tt, fo, sep = "\t", row.names = F, col.names = T, quote = F)

# Phenotype
phs = unique(ti$Phenotype)

fp = file.path(dirw, "11.phenotype.tsv")
tp = read.table(fp, header = T, sep = "\t", as.is = T)
sum(phs %in% tp$sname)

# Strain
## process entry data
fi = file.path(dirw, "URT_Entries04_05.RData")
attach(fi)
et = urt_entries

tt = data.frame()
for (yr in names(urt_entries)) {
	ds = urt_entries[[yr]]
	cols = colnames(ds)
	if("SeedSource" %in% cols) {seeds = ds$SeedSource} else {seeds = rep(NA, nrow(ds))}
	tts = data.frame(Year=rep(yr, nrow(ds)), Ent=ds$Ent, Strain=ds$Strain, Parentage=ds$Parentage, SeedSource=seeds, GenComp=ds$GenComp, UniqueTraits=ds$UniqueTraits, Trial=ds$Trial, PreviousTesting=ds$PreviousTesting, stringsAsFactors=F)
	tt = rbind(tt, tts)
}

tt$Year = substr(tt$Year, 2, 5)
tt$Parentage[tt$Parentage=='na'] = ''
tt[is.na(tt)] = ''
grp = dplyr::group_by(tt, Strain)
tt2 = dplyr::summarise(grp, 
	Year=paste(unique(Year[Year!='']), collapse=", "), 
	Parentage=Parentage[1], nParentage=length(unique(Parentage)), 
	SeedSource=paste(unique(SeedSource[SeedSource!='']), collapse=", "), 
	GenComp=paste(unique(GenComp[GenComp!='']), collapse=", "), 
	UniqueTraits=paste(unique(UniqueTraits[UniqueTraits!='']), collapse=", "), 
	Trial=paste(unique(Trial[Trial!='']), collapse=", "), 
	PreviousTesting=paste(unique(PreviousTesting[PreviousTesting!='']), collapse=", "))
table(tt2$nParentage)
tt = tt2[,-4]
#pnt = tt$Parentage
#idxs = grep("[^([{] ?[xX] ?[^([]", pnt)

## process all_ped data
fi = file.path(dirw, "all_ped_04_15.RData")
attach(fi)
tp = unique(all_ped[,c(-2,-3)])
tp[is.na(tp)] = ''
grp = dplyr::group_by(tp, Strain)
tp = dplyr::summarise(grp, Female=Female[1], Male=Male[1], 
	Synonyms=paste(unique(Synonyms[Synonyms!='']), collapse=", "), 
	Comments=paste(unique(Comments[Comments!='']), collapse=", "))

tm = merge(tp, tt, by='Strain', all.x=T)
tm[is.na(tm)]=''

strains = unique(ti$Strain)
length(strains)
strains_unk = strains[! strains %in% tm$Strain]
length(strains_unk)

pnt = unique(c(tm$Male, tm$Female))
pnt = pnt[pnt!='']
pnt_unk = pnt[! pnt %in% tm$Strain]
length(pnt)
length(pnt_unk)

to1 = data.frame(Strain=unique(c(strains_unk, pnt_unk)), Female='', Male='', Synonyms='', Comments='', Year='', Parentage='', SeedSource='', GenComp='', UniqueTraits='', Trial='', PreviousTesting='')
to = rbind(tm, to1)
to = cbind(id = 1:nrow(to), to)

tos = to[,1:2]; colnames(tos) = c('pid','pstrain')
tf1 = merge(to, tos, by.x='Female', by.y='pstrain', all.x=T)
colnames(tf1)[ncol(tf1)] = 'maternal_id'
tf2 = merge(tf1, tos, by.x='Male', by.y='pstrain', all.x=T)
colnames(tf2)[ncol(tf2)] = 'paternal_id'

to = tf2[,c(3,4,1,2,5:15)]
colnames(to) = c('id','name','maternal','paternal','alias_pi','alias_experiment',
	'year','parentage','seedsource','generation','unique_traits','trial',
	'previous_testing','maternal_id','paternal_id')
to = to[order(to$id),]
fo = file.path(dirw, "11.strain.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = T, na='')

##### process real data
fi = file.path(dirw, "10.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
colnames(ti) = c("year","test","phenotype","strain","location","value")
fl = file.path(dirw, "11.location.tsv")
tl = read.table(fl, header = T, sep = "\t", as.is = T)
ft = file.path(dirw, "11.test.tsv")
tt = read.table(ft, header = T, sep = "\t", as.is = T)
fr = file.path(dirw, "11.strain.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
fp = file.path(dirw, "11.phenotype.tsv")
tp = read.table(fp, header = T, sep = "\t", as.is = T)

## test_location
ttl = unique(ti[,c("year", "test", "location")])
ttl2 = merge(ttl, tl[,1:2], by.x='location',by.y='name')
colnames(ttl2)[ncol(ttl2)] = 'location_id'
stopifnot(nrow(ttl) == nrow(ttl2))
ttl3 = merge(ttl2, tt, by.x=c("year","test"), by.y=c("year","name"))
colnames(ttl3)[ncol(ttl3)] = 'test_id'
stopifnot(nrow(ttl) == nrow(ttl3))

ttl = cbind(id=1:nrow(ttl3), ttl3[order(ttl3$test_id, ttl3$location_id),])

to = cbind(ttl[,c(1,6,5)], row_spacing='', rows_per_plot='', yield_cv='', yield_lsd='', note='')
colnames(to)[2:3] = c("test", "location")
fo = file.path(dirw, "15.test_location.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = T, na='')

## trial (test_location_strain)
ta = unique(ti[,c("year", "test", "location", "strain")])
ta2 = merge(ta, ttl[,1:4], by=c("year","test","location"))
colnames(ta2)[ncol(ta2)] = 'test_location_id'
stopifnot(nrow(ta) == nrow(ta2))
ta3 = merge(ta2, tr[,1:2], by.x='strain', by.y='name')
colnames(ta3)[ncol(ta3)] = 'strain_id'
stopifnot(nrow(ta) == nrow(ta3))

ta = cbind(id=1:nrow(ta3), ta3[order(ta3$test_location_id, ta3$strain_id),])

to = cbind(ta[,c(1,6,7)], planting_date='', maturity_date='')
colnames(to)[2:3] = c("test_location", "strain")
fo = file.path(dirw, "16.trial.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = T, na='')

## observation
#ti = ti[!duplicated(ti[,1:5]),]
tv = unique(ti[,c("year", "test", "location", "strain", "phenotype")])
stopifnot(nrow(tv)==nrow(ti))

tv2 = merge(ti, ta[,1:5], by=c("year","test","location","strain"))
colnames(tv2)[ncol(tv2)] = 'trial_id'
stopifnot(nrow(ti) == nrow(tv2))
tv3 = merge(tv2, tp[,1:2], by.x='phenotype', by.y='name')
colnames(tv3)[ncol(tv3)] = 'phenotype_id'
stopifnot(nrow(ti) == nrow(tv3))

tn = cbind(id=1:nrow(tv3), tv3[order(tv3$trial_id, tv3$phenotype_id),])

to = tn[,c(1,8,9,7)]
colnames(to)[2:3] = c("trial","phenotype")
fo = file.path(dirw, "17.observation.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = T, na='')


         2YearMean          3YearMean          4YearMean        Descriptive
               188                 83                 14                  3
DescriptiveDisease            Entries             Height            Lodging
               408                220                414                414
          Maturity                Oil            Protein    RegionalSummary
               414                407                407                414
       SeedQuality           SeedSize           YieldBuA          YieldRank
               372                414                414                414
               
             BB              BSR        BSR%Incid         BSRPlant
               6                6               13               51
         BSR%Sev          BSRStem             BTSa        Chlorosis
              13               51               11              414
       Emergence               FE            FEarx             FELS
              57                2               82                1
       GreenStem         HardSeed           Mottle           NSC010
             258               63                1                6
        PhytoRot         PhytoTol               PM       PRPhytoTol
               2               23                1                1
         PRRace1          PRRace4          PRRace7               PS
              42              199              274                1
             PSa              PSB             P&SB          RootRot
             192               94               71               19
   RootRotRace25              SCL              SDS            SDSDI
               9                9              129               19
           SDSDS            SDSI%            SDSR6        SDSR6Date
              18                5                9                4
         SDSRank         SDSRDate             SDSS          SDSTest
               3                4                5                4
        SeedGerm       Shattering  Shattering10Oct Shattering30Sept
              28              392                1                1
  Shattering3Oct   Shattering9Oct              SMV            Stand
               1                1               32                7
      StemCanker
               1

require(dplyr)
require(ggmap)
require(xlsx)

dirw = file.path(Sys.getenv("misc2"), "urt")

### process raw data
fi = file.path(dirw, "urt_data.RData")
x = load(fi)
#fi = file.path(dirw, "NotesSummariesFromData.RData")
#y = load(fi)

### extract phenotypic data
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
		#if("Mean" %in% colnames(tp1)) { tp1 = tp1[,-which(colnames(tp1)=="Mean")] }
		if(ph == "Maturity") { for (k in 2:ncol(tp1)) { tp1[,k] = as.character(tp1[,k]) }}
		if(colnames(tp1)[1] != 'Strain') {
			cat(tnames[i], "not start with Strain: ", colnames(tp1)[1], "\n")
			colnames(tp1)[1] = "Strain"
		}
		tp2 = reshape(tp1, direction='long', varying=list(2:ncol(tp1)), idvar='Strain', timevar="Location", v.names='Value', times=colnames(tp1)[2:ncol(tp1)])
		tp3 = cbind.data.frame(Year=yr, Test=ex, Phenotype=ph, tp2, stringsAsFactors=F)
		tp3$Value = as.character(tp3$Value)
		tp = rbind(tp, tp3)
	} else if(ph %in% phs2) {
	} else if(ph %in% phs3) {
		tp1 = ud[[i]]
		if(colnames(tp1)[1] != 'Strain') {
			cat(tnames[i], "not start with Strain: ", colnames(tp1)[1], "\n")
			colnames(tp1)[1] = "Strain"
		}
		if(colnames(tp1)[2] != 'DescriptiveCode') {
			cat(tnames[i], "not start with DescriptiveCode: ", colnames(tp1)[2], "\n")
		} else {
			ph = 'DescriptiveCode'
			tp3 = data.frame(Year=yr, Test=ex, Phenotype=ph, Strain = tp1$Strain, Location='Mean', Value=as.character(tp1[,2]), stringsAsFactors = F)
			tp = rbind(tp, tp3)
		}
		if(ncol(tp1) <= 2) next
		zz = strsplit(colnames(tp1)[3:ncol(tp1)], split="[_]")
		for (j in 1:length(zz)) {
			ary = zz[[j]]
			if(length(ary) < 2) next
			if(length(ary) > 2) {
				cat(tnames[i], " col error: ", colnames(tp1)[j+2], "\n")
				if(ary[2] == '' | ary[2] == ary[3]) {
					ary[2] = ary[3]
				} else {
					stop("fatal error: unknown col")
				}
			}
			ph = ary[1]; lc = ary[2]
			tp3 = data.frame(Year=yr, Test=ex, Phenotype=ph, Strain = tp1$Strain, Location=lc, Value=as.character(tp1[,2+j]), stringsAsFactors = F)
			tp = rbind(tp, tp3)
		}
	} else {
		cat("unknown phenotype: ", ph, "\n")
	}
}

fo = file.path(dirw, "01.tsv")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### extract trial meta-data
pl = planting_maturity; yp = yield_plots

tpl = data.frame()
for (i in 1:length(pl)) {
	zz = strsplit(names(pl)[i], split="[_]")[[1]]
	yr = zz[1]; ex = zz[2]
	t1 = pl[[i]]
	for (j in 1:ncol(t1)) { t1[,j] = as.character(t1[,j]) }
	t2 = t1[(!is.na(t1[,1]) & t1[,1] != '') | (!is.na(t1[,2]) & t1[,2] != 0),]
	t3 = cbind.data.frame(Year=yr, Test=ex, Location=rownames(t2), t2, stringsAsFactors = F)
	t4 = reshape(t3, direction='long', varying=list(4:ncol(t3)), idvar=c("Year", "Test", "Location"), timevar="Meta", v.names='Value', times=colnames(t3)[4:ncol(t3)])
	tpl = rbind(tpl, t4)
}

typ = data.frame()
for (i in 1:length(yp)) {
	zz = strsplit(names(pl)[i], split="[_]")[[1]]
	yr = zz[1]; ex = zz[2]
	t1 = yp[[i]]
	colnames(t1)[1] = "Meta"
	t2 = t1#[t1$Meta != 'LocationMean',]
	for (j in 2:ncol(t2)) { t2[,j] = as.character(t2[,j]) }
	t3 = cbind.data.frame(Year=yr, Test=ex, t2, stringsAsFactors = F)
	t4 = reshape(t3, direction='long', varying=list(4:ncol(t3)), idvar=c("Year", "Test", "Meta"), timevar="Location", v.names='Value', times=colnames(t3)[4:ncol(t3)])
	typ = rbind(typ, t4)
}
typ = typ[,c(1,2,4,3,5)]
unique(typ$Meta)

tt1 = rbind(tpl, typ)
tt2 = reshape(tt1, direction = 'wide', idvar=c("Year","Test","Location"), timevar="Meta")
colnames(tt2)[4:ncol(tt2)] = gsub("Value.", "", colnames(tt2)[4:ncol(tt2)])

fo = file.path(dirw, "02.trial.meta.tsv")
write.table(tt2, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

##### data check
fi = file.path(dirw, "01.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
fr = file.path(dirw, "02.trial.meta.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)

# check year
ti$Year = as.numeric(substr(ti$Year, 2, 5))
tr$Year = as.numeric(substr(tr$Year, 2, 5))
unique(ti$Year)
unique(tr$Year)
identical(unique(ti$Year), unique(tr$Year))

# check phenotypes
phs = unique(ti$Phenotype)
phs

ph_map = c(
		'Shattering30Sept' = 'Shattering', 
		'Shattering9Oct'   = 'Shattering',
		'Shattering3Oct'   = 'Shattering',
		'Shattering10Oct'  = 'Shattering'
)
idxs = which(ti$Phenotype %in% names(ph_map))
ti$Phenotype[idxs] = ph_map[ti$Phenotype[idxs]]

phs = unique(ti$Phenotype)
phs
#write(phs, file.path(dirw, "03.phenotype.txt"))

# check locations
lcs = unique(ti$Location)
lcs

lc_map = c(
	'UllinILDS' = 'UllinIL',
	'Eldorado??' = 'EldoradoIL',
	'NWBranch' = 'NorthwestBranchMD',
	'FayettevilleirrigatedAR' = 'FayettevilleAR',
	'PinetreeirrigatedAR' = 'PineTreeAR',
	'LafayetteINLafayetteIN' = 'LafayetteIN',
	'PinetreeirrigatedAR' = 'PineTreeAR'
)
idxs = which(ti$Location %in% names(lc_map))
ti$Location[idxs] = lc_map[ti$Location[idxs]]
lcs = unique(ti$Location)

idxs = which(tr$Location %in% names(lc_map))
tr$Location[idxs] = lc_map[tr$Location[idxs]]
lcs2 = unique(tr$Location)
lcs_more = lcs2[! lcs2 %in% lcs]
if(length(lcs_more) > 0) lcs = c(lcs, lcs_more)

pre = sapply(lcs, jj <- function(x) substr(x, start=1, stop=nchar(x)-2))
suf = sapply(lcs, jj <- function(x) substr(x, start=nchar(x)-1, stop=nchar(x)))
idxs = which(suf %in% state.abb)
ts1 = data.frame(name = lcs[idxs], city = pre[idxs], state = suf[idxs], stringsAsFactors = F)

idxs = which(! suf %in% state.abb)
lcsc = lcs[idxs]
pre = sapply(lcsc, jj <- function(x) substr(x, start=1, stop=nchar(x)-3))
suf = sapply(lcsc, jj <- function(x) substr(x, start=nchar(x)-2, stop=nchar(x)))
ts2 = data.frame(name = lcsc, city = pre, state = suf, stringsAsFactors = F)
ts2$city[ts2$name == 'Mean'] = "Mean"; ts2$state[ts2$name == 'Mean'] = "Mean"

tl2 = rbind(ts1, ts2)
tl2 = tl[order(tl$state, tl$city),]

fl = file.path(dirw, "04.location.tsv")
tl = read.table(fl, sep = "\t", as.is = T, header = T)
stopifnot(tl2$name %in% tl$name)
#get_geocode(tl2, fl)
locmap = tl$location; names(locmap) = tl$name

get_geocode <- function(tl2, fo) {
	states=tl2$state; cities=tl2$city
	states[states=='QUE'] = 'QUEBEC'
	lonlat = geocode(sprintf("%s %s", cities, states))
	tl = cbind.data.frame(tl2, location = sprintf("%s_%s", tl2$city, tl2$state), longitude=lonlat[,1], latitude=lonlat[,2])
	write.table(tl, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

##
stopifnot(ti$Location %in% tl$name)
ti$Location = locmap[ti$Location]
fo = file.path(dirw, "10.tsv")
write.table(ti, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##
stopifnot(tr$Location %in% tl$name)
tr$Location = locmap[tr$Location]
tr = cbind.data.frame(tr, Trial = sprintf("%s_%s_%s", tr$Test, tr$Year, tr$Location), stringsAsFactors = F)
tr = tr[,c(1:3,12,4:11)]

tit = unique(ti[,c('Year','Test','Location')])
tit = cbind.data.frame(tit, Trial = sprintf("%s_%s_%s", tit$Test, tit$Year, tit$Location), stringsAsFactors = F)
tit2 = tit[!tit$Trial %in% tr$Trial,]
tre = cbind.data.frame(tit2, DatePlanted = sprintf("%d-01-01", tit2$Year), DaystoMature = 1, C.V. = NA, L.S.D. = NA, RowSp. = NA, Rows.Plot = NA, Reps = NA, LocationMean = NA, stringsAsFactors = F)

tr2 = rbind(tr, tre)
tr2 = tr2[order(tr2$Year, tr2$Test, tr2$Location),]
fo = file.path(dirw, "11.trial.meta.tsv")
write.table(tr2, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### write trial table for upload
fr = file.path(dirw, "11.trial.meta.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
fl = file.path(dirw, "04.location.tsv")
tl = read.table(fl, header = T, sep = "\t", as.is = T)

tr2 = merge(tr, tl[,c('location','longitude','latitude')], by.x = 'Location', by.y = 'location')
stopifnot(nrow(tr) == nrow(tr2))

desc = sprintf("RowSp.: %d, Rows/Plot: %d", tr2$RowSp., tr2$Rows.Plot)
pdates = as.Date(tr2$DatePlanted, "%Y-%m-%d")
hdates = pdates + tr2$DaystoMature
pdates = format(pdates, "%m/%d/%Y")
hdates = format(hdates, "%m/%d/%Y")
pdates[is.na(pdates)] = sprintf('1/1/%d', tr2$Year[is.na(pdates)])
hdates[is.na(hdates)] = sprintf('1/1/%d', tr2$Year[is.na(hdates)])
pdates = sub("0?(.+)/0?(.+)/(.+)", "\\1/\\2/\\3", pdates)
hdates = sub("0?(.+)/0?(.+)/(.+)", "\\1/\\2/\\3", hdates)
tr8 = cbind.data.frame(tr2[,c(4,2:3,1,13:14)], collab = 'Aaron Lorenz', trialdesc = '', pdate = pdates, hdate = hdates, weather = '', green = 'no', rate = '', design = '', nentry = '', nrep = '', plotsize = '', harea = '', irri = 'no', other = desc)
tr9 = t(tr8)
colnames(tr9) = sprintf("T%04d", 1:ncol(tr9))
rownames(tr9) = c(
"Trial Name", 
"Trial Year", 
"Experiment Name", 
"Location", 
"Latitude of field", 
"Longitude of field", 
"Collaborator", 
"Trial description", 
"Planting date", 
"Harvest date", 
"Begin weather date", 
"Greenhouse trial? (yes or no)", 
"Seeding rate (seeds/m2)", 
"Experimental design", 
"Number of entries", 
"Number of replications", 
"Plot size (m2)", 
"Harvested area (m2)", 
"Irrigation (yes or no)", 
"Other remarks")
tr9[,1:2]

psize = 2200
npart = ceiling(ncol(tr9)/psize)
for (i in 1:npart) {
fo = sprintf("%s/11.trial.meta/trial.%d.xlsx", dirw, i)
to = tr9[,((i-1)*psize+1):min(i*psize, ncol(tr9))]
wb <- createWorkbook()
sheet1 <- createSheet(wb, "Sheet1")
rows   <- createRow(sheet1, 1:4)
cells  <- createCell(rows, colIndex=1:2)
setCellValue(cells[[1,1]], "Trial Submission Form")
setCellValue(cells[[2,1]], "Template version")
setCellValue(cells[[2,2]], "20Dec2016")
setCellValue(cells[[3,1]], "Crop")
setCellValue(cells[[3,2]], "soybean")
setCellValue(cells[[4,1]], "Breeding Program Code")
setCellValue(cells[[4,2]], "URT")
addDataFrame(to, sheet1, col.names = T, row.names = T, startRow = 5, startColumn = 1)
#write.xlsx(x = tr9[,1:20], file = fo, sheetName = "Sheet1", row.names = T, col.names = T)
saveWorkbook(wb, fo)
}


# Strain
fi = file.path(dirw, "10.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
strains = unique(ti$Strain)
## process entry data
fi = file.path(dirw, "URT_Entries04_05.RData")
x2 = load(fi)
et = urt_entries


tnames = names(ud)
yrs = sapply(strsplit(tnames, split="[_]"), "[", 1)
exs = sapply(strsplit(tnames, split="[_]"), "[", 2)
phs = sapply(strsplit(tnames, split="[_]"), "[", 3)
idxs = which(phs == "Entries")
tt = data.frame()
for (idx in idxs) {
	ds = ud[[idx]]; yr = yrs[idx]; test = exs[idx]
	cols = colnames(ds)
	n = nrow(ds)
	if("SeedSource" %in% cols) {seeds = ds$SeedSource} else {seeds = rep(NA, n)}
	if("PreviousTesting" %in% cols) {ptests = ds$PreviousTesting} else {ptests = rep(NA, n)}
	if("UniqueTraits" %in% cols) {utrs = ds$UniqueTraits} else {utrs = rep(NA, n)}
	tts = data.frame(Year=rep(yr, nrow(ds)), Ent=ds$Ent, Strain=ds$Strain, Parentage=ds$Parentage, SeedSource=seeds, GenComp=ds$GenComp, UniqueTraits=utrs, Trial=rep(test, nrow(ds)), PreviousTesting=ptests, stringsAsFactors=F)
	tt = rbind(tt, tts)
}

tt$Year = substr(tt$Year, 2, 5)
tt$Parentage[tt$Parentage=='na'] = ''
tt[is.na(tt)] = ''

length(strains)
strainsn = unique(ti$Strain[ti$Year >= 2004 & ti$Year <= 2015])
strainsu = strainsn[!strainsn %in% tt$Strain]
length(strainsu)
fo = file.path(dirw, "21.parents.txt")
write(unique(tt$Parentage), fo)
fo = file.path(dirw, "21.entries.tsv")
write.table(tt, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = "")

grp = dplyr::group_by(tt, Strain)
tt = tt2[,-4]
tt2 = dplyr::summarise(grp, 
	Year=Year[1], Parentage=Parentage[1], SeedSource=SeedSource[1], 
	GenComp=GenComp[1], UniqueTraits=UniqueTraits[1], Trial=Trial[1], 
	PreviousTesting=PreviousTesting[1])
ti2 = ti[ti$Year >= 2004 & ti$Year <= 2015,]
sum(! unique(ti2$Strain) %in% tt2$Strain)

tx = data.frame(
	name = tt2$Strain, program = "URT", aliases = "", acc = "", 	
	pedigree = tt2$Parentage, generation = tt2$GenComp, species = "G.max", 
	comments = tt2$UniqueTraits, stringsAsFactors = F)
genmap = 1:9
names(genmap) = sprintf("F%d", 1:9)
tx$generation[!tx$generation %in% names(genmap)] = 'F1'
tx$generation = genmap[tx$generation]


f25 = sprintf("%s/25.strain.xlsx", dirw)
wb <- createWorkbook()
sheet1 <- createSheet(wb, "Sheet1")
rows   <- createRow(sheet1, 1:5)
cells  <- createCell(rows, colIndex=1:8)
setCellValue(cells[[1,1]], "Line Submission Form")
setCellValue(cells[[2,1]], "Version")
setCellValue(cells[[2,2]], "30Sep16")
setCellValue(cells[[3,1]], "*Crop")
setCellValue(cells[[3,2]], "Soybean")
setCellValue(cells[[4,1]], "*Line Name")
setCellValue(cells[[4,2]], "*Breeding Program")
setCellValue(cells[[4,3]], "Aliases")
setCellValue(cells[[4,4]], "GRIN Accession")
setCellValue(cells[[4,5]], "Pedigree")
setCellValue(cells[[4,6]], "*Filial Generation")
setCellValue(cells[[4,7]], "*Species")
setCellValue(cells[[4,8]], "Comments")
setCellValue(cells[[5,3]], "comma separated values")
setCellValue(cells[[5,5]], "Purdy notation")
setCellValue(cells[[5,6]], "0 - 9")
setCellValue(cells[[5,7]], "G.max / G.soja")
addDataFrame(tx, sheet1, col.names = F, row.names = F, startRow = 6, startColumn = 1, showNA = F)
saveWorkbook(wb, f25)


## write phenotype table for upload
options(java.parameters = "-Xmx10000m")
fi = file.path(dirw, "10.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tx = cbind.data.frame(ti, Trial = paste(ti$Test, ti$Year, ti$Location, sep = "_"), stringsAsFactors = F)[,c(7,4,3,6)]

fr = file.path(dirw, "11.trial.meta.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
tr2 = tr[!is.na(tr$LocationMean) | !is.na(tr$C.V.), c(4,12,7:8,11)]
tr2 = cbind(tr2, ste = (tr2$L.S.D. / 1.96) ^ 2 / 2)
tr2 = tr2[,c(1,2,6,5)] 
colnames(tr2)[2:ncol(tr2)] = c("*Trial Mean", "*Std. Error", "*Replications")
tr3 = reshape(tr2, direction = "long", idvar = c("Trial"), varying = list(2:ncol(tr2)), timevar = 'Strain', times = colnames(tr2)[2:ncol(tr2)], v.names = 'Value')
tr3 = tr3[order(tr3$Trial),]
tr4 = cbind.data.frame(tr3, Phenotype = "YieldBuA", stringsAsFactors = F)[,c(1,2,4,3)]

to = rbind(tx, tr4)
to = tx
yrs = sapply(strsplit(to$Trial, split = "_"), "[", 2)
for (yr in c(2004:2006)) {
	tos = to[yrs == yr,]
#	tos$Value[tos$Value == ''] = "NULL"
	t31 = reshape(tos, direction = 'wide', timevar = "Phenotype", idvar = c("Trial", "Strain"))
	colnames(t31)[-c(1:2)] = sapply(strsplit(colnames(t31)[-c(1:2)], split = "[.]"), "[", 2)
	t31 = cbind.data.frame(t31[1:2], check = 0, t31[-c(1:2)], stringsAsFactors = F)
	t31 = t31[order(t31$Trial, t31$Strain, method = 'radix'),]
	#write.table(t32, f32, sep = "\t", row.names = F, col.names = T, quote = F, na='')
	f31 = sprintf("%s/31.pheno/%s.xlsx", dirw, yr)

colnames(t31)[1:3] = c("*Trial Code", "*Line Name", "*Check")
wb <- createWorkbook()
sheet1 <- createSheet(wb, "Sheet1")
rows   <- createRow(sheet1, 1:2)
cells  <- createCell(rows, colIndex=1:2)
setCellValue(cells[[1,1]], "Phenotype Results")
setCellValue(cells[[2,1]], "Crop")
setCellValue(cells[[2,2]], "soybean")
addDataFrame(t31, sheet1, col.names = T, row.names = F, startRow = 3, startColumn = 1, showNA = T, characterNA = "NULL")
saveWorkbook(wb, f31)
}


#------obsolete
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
#pnt = tt$Parentage
#idxs = grep("[^([{] ?[xX] ?[^([]", pnt)

## process all_ped data
fi = file.path(dirw, "all_ped_04_15.RData")
x3 = load(fi)
tp = unique(all_ped[,c(-2,-3)])
tp[is.na(tp)] = ''
grp = dplyr::group_by(tp, Strain)
tp = dplyr::summarise(grp, Female=Female[1], Male=Male[1], 
	Synonyms=paste(unique(Synonyms[Synonyms!='']), collapse=", "), 
	Comments=paste(unique(Comments[Comments!='']), collapse=", "))

length(strains)
strainsu = strains[!strains %in% tt$Strain]
length(strainsu)
strainsu = strains[!strains %in% tt$Strain & !strains %in% ap$Strain]
length(strainsu)
head(strainsu, 50)


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

tx = data.frame(name=to$name, program="URT", aliases=to$alias_pi, acc="", pedigree=to$parentage, generation=to$generation, species="Y", comments=to$unique_traits, stringsAsFactors = F)
genmap = 1:9
names(genmap) = sprintf("F%d", 1:9)
tx$generation[!tx$generation %in% names(genmap)] = 'F1'
tx$generation = genmap[tx$generation]
fo = file.path(dirw, "11.strain.t3.tsv")
write.table(tx, fo, sep = "\t", row.names = F, col.names = T, quote = F, na='')


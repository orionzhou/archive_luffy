

dirw = file.path(Sys.getenv("misc2"), "urt")

fi = file.path(dirw, "urt_data.RData")

attach(fi)
ud = urt_data

tnames = names(ud)
yrs = sapply(strsplit(tnames, split="[.]"), "[", 1)
exs = sapply(strsplit(tnames, split="[.]"), "[", 2)
exs = sapply(strsplit(exs, split="[_]+"), "[", 1)
phs = sapply(strsplit(tnames, split="[_]+"), "[", 2)

cnames = unlist(sapply(ud, colnames))
res = strsplit(cnames, "[_]+")
idxs_di = which(sapply(res, length) > 1)
dnames = sapply(res[idxs_di], "[", 1)
lnames = cnames
lnames[idxs_di] = sapply(res[idxs_di], "[", 2)

         2YearMean          3YearMean          4YearMean        Descriptive
               188                 80                 14                  4
DescriptiveDisease            Entries             Height            Lodging
               405                103                411                411
          Maturity                Oil            Protein    RegionalSummary
               411                404                404                411
       SeedQuality           SeedSize           YieldBuA          YieldRank
               369                411                411                411

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

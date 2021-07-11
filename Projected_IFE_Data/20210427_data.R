setwd("~/Dropbox/InteractivefixedEffect/2020_Projection_IFE/Code/Application")

#Read in data from Penn World Table
df1 = read.csv2("pwt100.csv")
df1 = df1[c("countrycode","country","year","rgdpe","rgdpo","pop",
       "csh_c","csh_g","csh_i","pl_i")]
df1 = subset(df1,year>1969)
df1 = na.omit(df1)

#Read in WDI data
df2 = read.csv("wdi2.csv")[,-2]
names(df2) = c("year","country","countrycode","young","fert","life")
df2 = transform(df2,young=as.numeric(as.character(young)))
df2 = transform(df2,fert=as.numeric(as.character(fert)))
df2 = transform(df2,life=as.numeric(as.character(life)))

#Merge both datasets
data = merge(df1,df2,by=c("country","countrycode","year"))
data = na.omit(data)
data = data[data$country %in% names(which(table(data$country)==50)), ]
data$gdp.capita = data$rgdpe/data$pop
#data$gdp.growth = c(0,diff(log(data$rgdpe/data$pop)))
#data$pop.growth = c(0,diff(log(data$pop)))

data$gdp.growth = c(0,(data[-1,]$gdp.capita-data[-dim(data)[1],]$gdp.capita)/data[-dim(data)[1],]$gdp.capita)
data$pop.growth = c(0,(data[-1,]$pop-data[-dim(data)[1],]$pop)/data[-dim(data)[1],]$pop)

data = subset(data,year>1970)
data = data[data$country %in% names(which(table(data$country)>=49)), ]

write.csv(data,file="data_merged.csv",row.names = F)

##### NEW ####
df1 = read.csv2("pwt100.csv")
df1 = df1[c("countrycode","country","year","rgdpe","rgdpo","pop",
            "csh_c","csh_g","csh_i","pl_i")]
df1 = subset(df1,year>1988)
df1 = na.omit(df1)
df2 = read.csv("wdi2.csv")[,-2]
names(df2) = c("year","country","countrycode","young","fert","life")
df2 = transform(df2,young=as.numeric(as.character(young)))
df2 = transform(df2,fert=as.numeric(as.character(fert)))
df2 = transform(df2,life=as.numeric(as.character(life)))
data = merge(df1,df2,by=c("country","countrycode","year"))
data = data[data$country %in% names(which(table(data$country)==31)), ]
data$gdp.capita = data$rgdpe/data$pop
data$gdp.growth = c(0,(data[-1,]$gdp.capita-data[-dim(data)[1],]$gdp.capita)/data[-dim(data)[1],]$gdp.capita)
data$pop.growth = c(0,(data[-1,]$pop-data[-dim(data)[1],]$pop)/data[-dim(data)[1],]$pop)
data = subset(data,year>1989)
data = na.omit(data)
data = data[data$country %in% names(which(table(data$country)>=30)), ]
write.csv(data,file="data_merged.csv",row.names = F)

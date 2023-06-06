library(rmapshaper)

smallerws <- ms_simplify(watersheds)

leaflet(smallerws) %>%
                            addPolygons(weight = 1) 
                                                                               smoothFactor = 0.5,
                                                                               fillOpacity = 0)

writeOGR(obj = smallerws, 
         layer = "smallerws", 
         dsn = "temp",
         driver = "ESRI Shapefile")

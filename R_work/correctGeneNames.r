correctGeneNames <- function(oldGeneNames){
	newGeneNames 	<- c(
		"alx4b", "dpf3", "csf1ra", "ednrba", "ets1", "fgfr3", "foxd3", "foxg1b", "foxo1a", "foxo1b", 
		"foxp4", "hbp1", "her9", "hmx1", "hmx4", "id2a", "impdh1b", "gapdh", "kita", "ltk", "mbpa", 
		"mc1r", "mitfa", "mlphb", "mycla", "myo5aa", "oca2", "otx2b", "pax3a", "pax3b", 
		"pax7a", "pax7b", "phox2bb", "pnp4a", "pmela", "rpl13", "slc24a5", "smad9", "snai1b", "sox10", "sox5", 
		"sox9b", "tfap2a", "tfap2e", "tfec", "tyr", "tyrp1b")
	corrGeneNames		<- oldGeneNames
	corrGeneNames[ which( corrGeneNames == "csf1r")] 	<- "csf1ra"
	corrGeneNames[ which( corrGeneNames == "ets1a")] 	<- "ets1"
	corrGeneNames[ which( corrGeneNames == "fgfr3_v2")] 	<- "fgfr3"
	corrGeneNames[ which( corrGeneNames == "mycl1a")] 	<- "mycla"
	corrGeneNames[ which( corrGeneNames == "otx2")] 	<- "otx2b"
	corrGeneNames[ which( corrGeneNames == "pax3_v2")] 	<- "pax3a"
	corrGeneNames[ which( corrGeneNames == "phox2b")] 	<- "phox2bb"
	corrGeneNames[ which( corrGeneNames == "silva")] 	<- "pmela"
	corrGeneNames[ which( corrGeneNames == "snail2")] 	<- "snai1b"
return( corrGeneNames)
}

# this function takes the string 'sourceStr' and returns the substring identified by regexpRes made by regexpr()

getRegExpSubStr <- function(sourceStr, regexprRes){		

ans <-  substr( sourceStr, regexprRes, regexprRes + attr(regexprRes, "match.length") -1)
ans
}


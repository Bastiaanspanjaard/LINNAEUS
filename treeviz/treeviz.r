do_install = function(){
	install.packages("treemap")
	install.packages("data.tree")
	install.packages("collapsibleTree")
}

rl = function(){
	library(data.tree)
	library(treemap)
	library(collapsibleTree)
}

#dev_tree = read.table('/Users/polivar/Downloads/bs/tree_B_dev_tree.txt', h=T)
#B_scartree = read.table('/Users/polivar/Downloads/bs/tree_B_scar_tree.txt', h=T)
#> head(B_scartree)
#  row.names V1    V2 Scar.acquisition
#1        26  7 37970               NA
#2        26  8 24841               NA
#3        26  1 26212            39244
#4        25  9 33957               NA
#5        27 10 24547               NA
#6        27 11 19955               NA

#hh = head(dev_tree)

generate_tree = function(df){
	rm(list=ls(envir=globalenv(), pattern='^n'), envir=globalenv())
	apply(df, 1, function(x){
		parent = paste0('nd', as.numeric(x[1]))
		child = paste0('nd', as.numeric(x[2]))
		if(length(x) ==4){
			scar_p = ifelse(is.na(x[3]), parent, x[3])
			scar_c = ifelse(is.na(x[4]), child, x[4])
		}
		if(length(x) ==3){
			scar_p = ' '
			scar_c = ifelse(is.na(x[3]), child, x[3])
		}
		if(!exists(parent)){ 
			eval_txt = sprintf('%s <<- Node$new("%s", name="%s")', parent, parent, scar_p)
	#		message(eval_txt)
			eval(parse(text=eval_txt))
		}
		if(!exists(child)){ 
			eval_txt = sprintf('%s <<- Node$new("%s", name="%s")', child, child, scar_c)
	#		message(eval_txt)
			eval(parse(text=eval_txt))
		}
		add_txt = sprintf('%s$AddChildNode(%s)', parent, child)
		eval(parse(text=add_txt))
	})
	return(eval(parse(text=sprintf('%s$root', ls(envir=globalenv(), pattern='^nd')[1]))))
}

#for( i in ls(pattern='^n')){eval(parse(text=sprintf('print(%s$level)',i)))}

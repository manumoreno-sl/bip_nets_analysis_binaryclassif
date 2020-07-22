# """
#         Written by Manuel Moreno , maamorenova@unal.edu.co, 2020.
#        
#         Dyads Variable Analysis for Categorical Variables Embedding        
# 		  R version 3.6
# """
## Dependencies ###########################################################
require(igraph)
require(magrittr)
require(dplyr)
# require(bipartite)
require(tnet)
require(ade4)
require(psych)
require(assortnet)

source('./r package/bivariate_moran.R')

## Functions for Creating Bipartite and their projections ###########################################################
norm_scaler=function(x){
  # q=quantile(x)[4]-quantile(x)[2]
  # x=x/q
  normscal=x/sum(x)
  return(normscal) #x
}

proj_wt_calc=function(x,cvar,type_num,dat){
	dat %>% 
  	dplyr::filter(!!sym(cvar) %in% x) %>% 
  	dplyr::summarize( #val1=x[1],val2=x[2],
  		t_wt=sum(t_wt),y_wt=sum(y_wt),
		n_ed=dplyr::n())
}

proj_wt_proc=function(proj,cvar_n,dat){
	# (g_uni$proj2,'2',dat)

	ed2=as.data.frame(do.call(rbind, strsplit(as_ids(E(proj)),"\\|")))
	e_attr=apply(ed2,MARGIN = 1,
				 FUN=proj_wt_calc,cvar=paste("cvar",cvar_n,sep=''),type_num=2,dat=dat)
	e_attrs=as.data.frame(do.call(rbind,e_attr))
	edge_attr(proj, "t_wt") <- e_attrs$t_wt
	edge_attr(proj, "y_wt") <- e_attrs$y_wt
	edge_attr(proj, "n_ed") <- e_attrs$n_ed
	return(proj)
}

proj_type=function(bip_net,method="overlap",bin=FALSE){
	# Method= overlap, jaccard, matching, pearson, yule
	# Bin= binarize? if false get raw metrics. If a numeric value, treshold to binarize
	inc_mat=as_incidence_matrix(bip_net)  # Extract the matrix

	if(method=="overlap"){
	#Overlap count##############################################################
 	projected_g <- bipartite_projection(g, multiplicity = TRUE)
 	g_proj1=g_proj1$proj1
 	g_proj1=g_proj2$proj2
 	} else if(method=="jaccard"){
	#Jaccard Similarity ########################################################
	proj1_jaccard <- as.matrix(dist.binary(inc_mat, method=1, upper=TRUE, diag = FALSE)) # Method #1 is "Jaccard Index"
	diag(proj1_jaccard)<-0
	proj2_jaccard <- as.matrix(dist.binary(t(inc_mat), method=1, upper=TRUE, diag = FALSE))
	diag(proj2_jaccard)<-0
	if(bin!=FALSE){
		proj1_jaccard <- ifelse(proj1_jaccard>bin, 1, 0)     # Binarize
		proj2_jaccard <- ifelse(proj2_jaccard>bin, 1, 0)     # Binarize
	} 
	g_proj1 <- graph_from_adjacency_matrix(proj1_jaccard,mode = "undirected")
	g_proj2 <- graph_from_adjacency_matrix(proj2_jaccard,mode = "undirected")

	} else if(method=="matching"){
	#Simple matching ###########################################################
	proj1_sm <- as.matrix(dist.binary(inc_mat, method=2, upper=TRUE, diag = FALSE))
	diag(proj1_sm)<-0
	proj2_jaccard <- as.matrix(dist.binary(t(inc_mat), method=2, upper=TRUE, diag = FALSE))
	diag(proj2_sm)<-0
	if(bin!=FALSE){
		proj1_sm <- ifelse(proj1_sm>bin, 1, 0)     # Binarize
		proj2_sm <- ifelse(proj2_sm>bin, 1, 0)     # Binarize
	} 
	g_proj1 <- graph_from_adjacency_matrix(proj1_sm,mode = "undirected")
	g_proj2 <- graph_from_adjacency_matrix(proj2_sm,mode = "undirected")

	} else if(method=="matching"){
	#Pearsons Correlation#######################################################
	proj2_corr <- as.matrix(cor(t(inc_mat)))
	diag(proj2_corr)<-0
	proj1_corr <- as.matrix(cor(inc_mat))
	diag(proj1_corr)<-0
	if(bin!=FALSE){
		proj1_corr <- ifelse(proj1_corr>bin, 1, 0)     # Binarize
		proj2_corr <- ifelse(proj2_corr>bin, 1, 0)     # Binarize
	} 
	g_proj1 <- graph_from_adjacency_matrix(proj1_corr,mode = "undirected")
	g_proj2 <- graph_from_adjacency_matrix(proj2_corr,mode = "undirected")

	} else if(method=="yule"){
	#Yule's Q###################################################################
	proj2_Q <-as.matrix(YuleCor(t(inc_mat))$rho)
	diag(proj2_Q)<-0
	proj1_Q <-as.matrix(YuleCor(inc_mat)$rho)
	diag(proj1_Q)<-0
	if(bin!=FALSE){
		proj1_Q <- ifelse(proj1_Q>bin, 1, 0)     # Binarize
		proj2_Q <- ifelse(proj2_Q>bin, 1, 0)     # Binarize
	}
	g_proj1 <- graph_from_adjacency_matrix(proj1_Q,mode = "undirected")
	g_proj2 <- graph_from_adjacency_matrix(proj2_Q,mode = "undirected")
	}
	return(list(proj1=g_proj1,proj2=g_proj2))
}

write_dat_binet=function(g_long,path,vars){
	name_vars=paste(paste(c("binet",vars),collapse = "_"),"dat",sep=".")
	write.table(g_long, file = paste(path,name_vars,sep="\\"), row.names=FALSE,col.names = FALSE,
		sep = "\t", quote = FALSE)
}

write_txt_binet=function(rosetta,path,vars){
	name_vars=paste(paste(c("binet",vars,"rosetta"),collapse = "_"),"dat",sep=".")
	write.table(rosetta, file = paste(path,name_vars,sep="\\"), row.names=FALSE, col.names = TRUE,
		sep = "\t", quote = FALSE)
}

bip_net=function(x,t_weight=NULL,y_weight,data,net_weight_type="y",calc_rank=TRUE,
				 damping=0.85,write_binet=FALSE,path=".\\"){
	print(x)
	# Create biparite
	# x= vector of length 2 with variables names of categorical columns
	# t_wt= name of time weight column
	# y_wt= name of weighted response variable column
	# data= dataset
	cvar1=x[1]
	cvar2=x[2]
	g_size=dim(data)[1]

	dat=claims %>% dplyr::select(all_of(c(cvar1,cvar2,t_weight,y_weight))) %>% 
		dplyr::rename(cvar1=all_of(cvar1),cvar2=all_of(cvar2),t_wt=all_of(t_weight),
			y_wt=all_of(y_weight)) 
	g <- graph_from_data_frame(
		dat %>%
		dplyr::group_by(cvar1,cvar2) %>% 
		#Setting Edge weights, related to t_wt
		dplyr::summarize(t_wt=sum(t_wt)/dim(.)[1],y_wt=sum(y_wt),n_ed=dplyr::n()), 
		directed = FALSE)
	V(g)$y_str=norm_scaler(igraph::strength(g,weights = E(g)$y_wt))
	V(g)$t_str=norm_scaler(igraph::strength(g,weights = E(g)$t_wt))

	if(calc_rank==TRUE){
	gotcha_rank=page_rank(g,directed = F,algo = "prpack",damping = damping,
		personalized = norm_scaler(V(g)$y_str*igraph::degree(g)),
		weights = E(g)$t_wt)$vector
	V(g)$grank=unname(gotcha_rank)
	}

  	if(bipartite.mapping(g)$res){
  		V(g)$type <- bipartite_mapping(g)$type
  		g_uni=bipartite_projection(g,remove.type=FALSE)
  		g_uni$proj1=proj_wt_proc(g_uni$proj1,'1',dat)
  		g_uni$proj2=proj_wt_proc(g_uni$proj2,'2',dat)
  	} else {
  		return("Non existing two-mode partition for these variables")
  	}

  	# if(net_weight_type=="y"){
  	# 	dat_bipk=dat %>% mutate(cg="i") %>% dplyr::select(cvar1,cvar2,cg,y_wt)	%>% rename(wt=y_wt)
  	# } else{
  	# 	dat_bipk=dat %>% mutate(cg="i") %>% dplyr::select(cvar1,cvar2,cg,t_wt)	%>% rename(wt=t_wt)
  	# }
  	# as_incidence_matrix(g)
  	# tweb=frame2webs(dat_bipk,varnames=c("cvar1","cvar2","cg","wt"),type.out="array")
  	# web=tweb[,,"i"]

  	# BINET requirement################
  	g_long=as_long_data_frame(g)
  	g_long_binet_form=g_long %>% 
  		mutate(from=paste("u",from,sep=""),to=paste("i",to,sep = "")) %>%
  		dplyr::select(from,to ,y_wt)
  	binet_rosetta=g_long %>% mutate(from=paste("u",from,sep="")) %>% dplyr::select(id=from,orig=from_name) %>% unique(.) %>%
  		bind_rows(g_long %>% mutate(to=paste("i",to,sep="")) %>%  dplyr::select(id=to,orig=to_name) %>% unique(.))
	
	if(write_binet==TRUE){
		print(c("writing BiNet to",path))
		write_dat_binet(g_long=g_long_binet_form,path=path,vars=x)
		write_txt_binet(rosetta=binet_rosetta,path=path,vars=x)
	}

  	if(bipartite.mapping(g)$res){
  		return(list(g_bip=g,g_uni=g_uni,g_bip_long_f=g_long_binet_form,vars=x,g_bip_names=binet_rosetta#,web=web
  			))
  	} else {
  		return(list(g_bip=g
  			# web=web
  			))
  	}
}


# Functions for measures  ###############################################################
# Global 

bip_density=function(g_bip){
	typeI_n=sum(V(g_bip)$type)
	####Networkx equivalent  https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/bipartite/basic.html#density
	n = gorder(g_bip)
    m = gsize(g_bip)
    nb = typeI_n#len(nodes)
    nt = n - nb
    if(m == 0){  # includes cases n==0 and n==1
        d = 0.0
    }else{
        # if B.is_directed():
        #     d = m / (2.0 * float(nb * nt))
        # else:
        d = m / (nb * nt)
     }
    return(d)
}

measures_glob=function(bip_net,z,alpha=0.05,vars,degree_moran){
	#bip_net= bip_net$g_bip
	#z- name of the weights attrubute
	A=bip_net[]

	###############################################################################################
	assort=assortment.continuous(A, vertex_attr(bip_net,z) , weighted = TRUE, SE = FALSE, M = 1)$r
	#Moran I####################################################################
	
	moran=moran_calc(A=A,vertex_attr(bip_net,z),autocorrelation=TRUE,alpha=alpha,degree=degree_moran)$glob # glob, loc, sig
	############################################################################
	ord = gorder(bip_net)
	siz = gsize(bip_net)
	nv = sum(V(bip_net)$type)#len(nodes)
	nu = ord - nv

	order=length(V(bip_net))
	# #- Size (number of edges Ne = |E|)
	size=length(E(bip_net))
	# #- Diameter (geodesic distance): 
	diameter=diameter(bip_net)
	# #- Average Distance: Average geodesic distance in the network
	avg_distance=average.path.length(bip_net) 
	# ## Measures of Interconnectedness
	# #- Density
	density=bip_density(bip_net)#edge_density(bip_net)
	# #- Maximum Matching
	max_bip_match_size=max_bipartite_match(bip_net)$matching_size
	max_bip_wmatch_size=max_bipartite_match(bip_net)$matching_weight
	# #- Spectral Bipartivity
	## In the works
	# #- Average Degree
	n_degree=igraph::degree(bip_net)
	u_deg_centrality=n_degree[V(bip_net)$type==FALSE]/nv
	v_deg_centrality=n_degree[V(bip_net)$type==TRUE]/nu
	avg_degree_var1=mean(u_deg_centrality)
	avg_degree_var2=mean(v_deg_centrality)
	
	# # Centralization
	
	centz_deg_u=centralize(u_deg_centrality,(nu-1)*(nu-2),TRUE)
	centz_deg_v=centralize(v_deg_centrality,(nv-1)*(nv-2),TRUE)
	# centr_bet=centr_betw(bip_net)$centralization
	# centr_clo=centr_clo(bip_net)$centralization
	
	bip_net_adj_w=as_adjacency_matrix(bip_net,attr = "t_wt",sparse = T)
	bip_adj=as.matrix(bip_net_adj_w[1:nu,(nu+1):(nu+nv)])
	W.svd <- svd(bip_adj,nu=1,nv=1)
	U_eig_score=abs(W.svd$u)
	V_eig_score=abs(W.svd$v)
	centz_eig_u=centralize(scores = U_eig_score,theoretical.max = nu-2,normalized = T)
	centz_eig_v=centralize(scores = V_eig_score,theoretical.max = nv-2,normalized = T)

	# # Standard Deviation
	# sd_deg=sd(igraph::degree(g))      # Business
	# sd_bet=sd(igraph::betweenness(g))
	# sd_clo=sd(igraph::closeness(g))
	# sd_eig=sd(igraph::evcent(g)$vector)
	return(data.frame(var1=vars[1],var2=vars[2],assort=assort,moranI=moran,order=order,size=size,diameter=diameter,
		avg_distance=avg_distance,density=density,max_bip_match_size=max_bip_match_size,
		max_bip_wmatch_size=max_bip_wmatch_size,centz_deg_u=centz_deg_u,centz_deg_v=centz_deg_v,avg_degree_var1=avg_degree_var1,avg_degree_var2=avg_degree_var2,
		centz_eig_u=centz_eig_u,centz_eig_v=centz_eig_v
		# centr_deg=centr_deg,centr_clo=centr_clo,centr_eig=centr_eig#,centr_bet=centr_bet
		))
}

#Node level 

bip_close=function(bip_net,nu,nv){
    n_closs=igraph::closeness(bip_net)
    u_closs0=n_closs[V(bip_net)$type==FALSE]
    v_closs0=n_closs[V(bip_net)$type==TRUE]
    u_closs=(nv+2*(nu-1))*u_closs0
    v_closs=(nu+2*(nv-1))*v_closs0
    return(list(u_closs=u_closs,v_closs=v_closs))
}

local_measures=function(bip_net,z,alpha=0.05,vars,degree_moran){
	ord = gorder(bip_net)
	siz = gsize(bip_net)
	nv = sum(V(bip_net)$type)#len(nodes)
	nu = ord - nv

	ranks=vertex_attr(bip_net,z)
	names=vertex_attr(bip_net,"name")
	################################################
	A=bip_net[]
	moran=moran_calc(A=A,ranks,autocorrelation=TRUE,alpha=alpha,degree=degree_moran)
	
	moran_l=moran$loc # 
	moran_s=moran$sig
	local_moran=data.frame(node=moran$v_name,moran_loc=moran_l,moran_sig=moran_s)
	#################################################
	bip_closeness=bip_close(bip_net,nu,nv)
	closeness_u=bip_closeness[[1]]
	closeness_v=bip_closeness[[2]]
	closeness=unname(c(closeness_u,closeness_v))
	local_closeness=data.frame(node=names,closeness=closeness)
	###################################################
	bip_net_adj_w=as_adjacency_matrix(bip_net,attr = "t_wt",sparse = T)
	bip_adj=as.matrix(bip_net_adj_w[1:nu,(nu+1):(nu+nv)])
	W.svd <- svd(bip_adj,nu=1,nv=1)
	U_eig_score=abs(W.svd$u)
	V_eig_score=abs(W.svd$v)
	eigen_cen=unname(c(U_eig_score,V_eig_score))
	local_eig=data.frame(node=names,eig_cen=eigen_cen)

	nodes_df=data.frame(var1=vars[1],var2=vars[2],node=names,z=ranks)	%>%
				left_join(local_moran,by=c("node")) %>%
				left_join(local_closeness,by=c("node")) %>%
				left_join(local_eig,by=c("node")) 

	return(nodes_df)#list(node=moran$v_name,local=moran_l,sig=moran_s))
}


# measures_cent_binary=function(){
# 	women_deg <- degree(jacc_women)
# 	women_bet <- betweenness(jacc_women)
# 	women_clos <- closeness(jacc_women)
# 	women_eig <- eigen_centrality(jacc_women)$vector
# 	women_cent_df <- data.frame(women_deg, women_bet, women_clos, women_eig)
# }

# measures_cent_weighted=function(){
# 	JW <- as.tnet(w_graph)
# 	Wdeg <- degree_w(JW)# [,2]
# 	Wbet <- betweenness_w(JW)# [,2]
# 	Wclos <- closeness_w(JW, gconly=FALSE) # [,2]
# 	W_cent <- data.frame(Wdeg, Wbet, Wclos)
# 	return(W_cent)
# }
# # Option 2: Analyze Each Mode Separately Using Two-Mode Metrics
# bip_degree=function(bip_net){
# 	NodeLabels <- V(bip_net)$name
# 	tm<-get.edgelist(bip_net, names=FALSE)
# 	mt <- tm[, c(2, 1)]
# 	deg_tm <- degree_tm(tm)
# 	deg_mt <- degree_tm(mt)
# 	return(list(deg_tm,deg_mt))
# }

# measures_net_loc=function(bip_net,mode="two"){
# 	if(mode=="two"){
# 		node_array=specieslevel(bip_net$web, level="both", index=c(
# 	  	# Not to include:'degree', 'interaction', 
# 	  	# Not useful for bipartite networks: 'betweenness',
# 	  	# Test and analyze
# 	  	"normalised degree","closeness", 'nestedrank',"PDI","species strength"
# 		), PDI.normalise=T)
# 		node_h=cbind(node=rownames(node_array$`higher level`),node_array$`higher level`)
# 		node_l=cbind(node=rownames(node_array$`lower level`),node_array$`lower level`)
# 		nodes=rbind(node_l,node_h) %>% dplyr::select(-closeness)
# 		return(nodes)

# 	} else if(mode=="one_l"){
# 		node_array=specieslevel(bip_net$web, level="lower", index=c(
# 	  	'degree', 'interaction', 
# 	  	'betweenness',
# 	  	"normalised degree","closeness", 'nestedrank',"PDI","species strength"
# 		), PDI.normalise=T)
# 		# node_h=cbind(node=rownames(node_array$`higher level`),node_array$`higher level`)
# 		node_l=cbind(node=rownames(node_array),node_array)
# 		# nodes=rbind(node_l,node_h) %>% dplyr::select(-closeness)
# 		return(node_l)
# 		#bipartite calcula el univariado, hacer formulas
# 	}
# }



# Functions for Plotting graphs ########################################
                       
get_loc=function(x){
  loc=x$loc_measures
  return(loc)
}

get_glob=function(x){
  glob=x$glob_measures
  return(glob)
}

bip_embed=function(x,data,t_w="time_weight",y_w='fraud_weight',calc_rank=TRUE, 
				   net_weight_type="t",damping=0.85,write_binet=FALSE,path=".\\",
				   alpha=0.05,degree_moran){

	g_list=bip_net(x=x,t_weight=t_w,y_weight=y_w,data=data,
		net_weight_type=net_weight_type,calc_rank=calc_rank,
		damping=damping,write_binet=write_binet,path=path)
	if(is.character(g_list)==TRUE){
		return("not a bipartie partition found")
	} else{
	glob_measures=measures_glob(bip_net=g_list$g_bip,"grank",alpha=alpha,g_list$vars,degree_moran)
	loc_measures=local_measures(bip_net=g_list$g_bip,"grank",alpha=alpha,g_list$vars,degree_moran)
	}
	# measures=c(global,local)
	return(list(g_list=g_list,glob_measures=glob_measures,loc_measures=loc_measures
		))
}

start_bip_build=function(nodes_vars,bip_embd_func,data,bine_folder,write_binet=FALSE,
	t_w="time_weight",y_w='fraud_weight',calc_rank=TRUE,damping=0.85,alpha=0.05,
	net_weight_type="t",degree_moran
	){
	require(tictoc)
	tic()
	comb_bip_bild=utils::combn(x=nodes_vars,m=2,simplify=F, FUN=bip_embd_func,data=data,write_binet=write_binet,path=bine_folder,
		t_w=t_w,y_w=y_w,calc_rank=calc_rank,damping=damping,alpha=alpha,
		net_weight_type=net_weight_type,degree_moran=degree_moran)
	tooc=toc()

	bip_guide=t(utils::combn(x=nodes_vars,m=2))
	colnames(bip_guide)=c("var1","var2")
	
	running_time=tooc$toc-tooc$tic
	
	comb_bip_bild_alt=comb_bip_bild[lapply(comb_bip_bild,is.character)!=TRUE]

	loc_df=bind_rows(lapply(comb_bip_bild_alt,get_loc))# .id = "column_label")
	glob_df=bind_rows(lapply(comb_bip_bild_alt,get_glob))# .id = "column_label")

	return(list(bip_list=comb_bip_bild,
			    all_glob_meas=glob_df,all_local_meas=loc_df,
			    bip_guide=bip_guide,time=running_time
	))
}


# Functions for Multipartite graphs ########################################


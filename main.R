setwd("/home/maddy/programs/blog post model")

# source("sims.R")
source("jacobian.R")
library(igraph)
library(network)

# jac<-getJacobian(steadystates[-1,2])
jac<-getJacobian(rep(1,27))
adjacencymatrix<-t(jac!=0)

net<-as.network(x=adjacencymatrix,directed=TRUE,
                matrix.type="adjacency",loops=TRUE)
net_noloops<-as.network(x=adjacencymatrix,directed=TRUE,
                        matrix.type="adjacency",loops=FALSE)

species<-c("IR","IRp","IRins","IRip","IRi","IRS1","IRS1p",
           "IRS1p307","IRS1307", "X","Xp","PKB","PKB308p","PKB473p",
           "PKB308p473p","mTORC1","mTORC1a","mTORC2","mTORC2a",
           "AS160","AS160p","GLUT4m","GLUT4","S6K","S6Kp","S6","S6p");

network.vertex.names(net)<-species
network.vertex.names(net_noloops)<-species


# indegree and outdegree are scientific names
indegree<-rowSums(adjacencymatrix)
outdegree<-colSums(adjacencymatrix)
total_deg<-indegree + outdegree

indegree_noloops<-rowSums(adjacencymatrix)-1
outdegree_noloops<-colSums(adjacencymatrix)-1
total_deg_noloops<-indegree_noloops + outdegree_noloops

jpeg("img/fullnet.jpg",width=1000,height=1000)
set.seed(1);
plot(net,displaylabels=FALSE,vertex.cex=total_deg/5,
     vertex.col=c("black","red")[1+(total_deg_noloops>6)],
     vertex.border=0,label.pos=1)
dev.off();


jpeg("img/net_noloops.jpg");
set.seed(1);
plot(net_noloops,displaylabels=TRUE,vertex.cex=total_deg_noloops/5,
     vertex.col=c("black","red")[1+(total_deg_noloops>6)],
     vertex.border=0,label.pos=1)
dev.off()

#igraph
ig<-graph_from_adjacency_matrix(adjacencymatrix,diag=FALSE,add.colnames = TRUE)

degree_df<-data.frame("vertices"=species,"deg_in"=degree(ig,mode="in"),"deg_out"=degree(ig,mode="out"),"deg_tot"=degree(ig))
degree_df <- degree_df %>% arrange(desc(deg_tot)) 

# degree distribution
ggplot()+geom_histogram(aes(degree_df$deg_tot),binwidth=1);

pdf_degree<-data.frame("degree"=0:10,"p"=degree_distribution(ig));
pdf_degree<-pdf_degree[c(-1,-2),];

jpeg("img/log_pdf.jpg");
ggplot()+geom_point(aes(x=log(pdf_degree$degree),y=log(pdf_degree$p)))+
  guides(size="none",linetype="none")+
  labs(x="",y="",title="Log transformed pdf")+
  theme_classic()+
  theme(text=element_text(size=20,face = "bold"))
dev.off();


bestfit1<-lm(log(pdf_degree$p) ~ log(pdf_degree$degree))

c=exp(bestfit1$coefficients[[1]]);
a=bestfit1$coefficients[[2]];
l=mean(degree_df$deg_tot); # Lamda in Poisson distribution
  
pdf_degree$powerlaw<-c*pdf_degree$degree^a;
pdf_degree$poisson<-exp(-l)*l^(pdf_degree$degree)/factorial(pdf_degree$degree);

jpeg("img/pdf_overlay.jpg");
gg<-ggplot(pdf_degree,aes(x=degree))
gg+geom_line(aes(y=p,color="Current System"),size=1.2)+
  geom_line(aes(y=powerlaw,color="Scale Free"),size=1.2)+
  geom_line(aes(y=poisson,color="Random network"),size=1.2)+
  guides(size="none",linetype="none")+
  labs(x="Degree")+
  ylim(0,1)+theme_classic()+
  theme(text=element_text(size=16,face="bold"),legend.box="vertical",legend.position="bottom");
dev.off();

loglikelihood_poisson=sum(log(pdf_degree$poisson[degree_df$deg_tot-1]));
loglikelihood_powerlaw=sum(log(pdf_degree$powerlaw[degree_df$deg_tot-1]));


# Ignoring the degree 2
bestfit2<-lm(log(pdf_degree$p[-1]) ~ log(pdf_degree$degree[-1]))

c=exp(bestfit2$coefficients[[1]]);
a=bestfit2$coefficients[[2]];

pdf_degree$powerlaw_nofirst<-c*pdf_degree$degree^a;
pdf_degree<-pdf_degree[-1,];

jpeg("img/pdf_overlay_noFirst.jpg");
gg<-ggplot(pdf_degree,aes(x=degree))
gg+geom_line(aes(y=p,color="Current System"),size=1.2)+
  geom_line(aes(y=powerlaw_nofirst,color="Scale Free"),size=1.2)+
  geom_line(aes(y=poisson,color="Random network"),size=1.2)+
  guides(size="none",linetype="none")+
  labs(x="Degree")+
  ylim(0,1)+theme_classic()+
  theme(text=element_text(size=16,face = "bold"),legend.position = "bottom");
dev.off();


loglikelihood_poisson_nooutlier=sum(log(pdf_degree$poisson[degree_df$deg_tot-2]));
loglikelihood_powerlaw_nooutlier=sum(log(pdf_degree$powerlaw[degree_df$deg_tot-2]));


library("deSolve")

SEIR.beta.gamma <- function(t, x, parms)
{
  S <- x[1] 
  E <- x[2]
  I <- x[3]
  R <- x[4]
  
  with(as.list(parms),
{ 
  dS <- mu - beta*S*I - mu*S 
  dE <- beta*S*I - sigma*E - mu*E
  dI <- sigma*E - (gamma+mu)*I
  dR <- (1-f)*gamma*I    # f = fatality rate
  out <- c(dS, dE, dI,dR) 
  list(out)
})
} 


SEIR.beta.gamma.N <- function(t, x, parms)
{
  # Same as SEIR.beta.gamma.N
  # but populations are in actual 
  # number of individuals, not proportions
  
  S <- x[1] 
  E <- x[2]
  I <- x[3]
  R <- x[4]
  with(as.list(parms),
{ 
  dS <- mu*N - beta*S*I/N - mu*S 
  dE <- beta*S*I/N - sigma*E - mu*E
  dI <- sigma*E - (gamma+mu)*I
  dR <- (1-f)*gamma*I    # f = fatality rate
  out <- c(dS, dE, dI,dR) 
  list(out)
})
} 


SIR.beta.gamma.N <- function(t, x, parms)
{
  
  S <- x[1] 
  I <- x[2]
  R <- x[3]
  
  with(as.list(parms),
{ 
  dS <- mu*N - beta*S*I/N - mu*S 
  dI <- beta*S*I/N - (gamma+mu)*I
  dR <- (1-f)*gamma*I    # f = fatality rate
  out <- c(dS, dI,dR) 
  list(out)
})
} 

control.transmission <- function(t, tau, k)
{
  # decrease (for transmission rate) with time
  # after time 'tau' and effort 'k'
  return(exp(-k*max(0,t-tau)))
}

control.transmission2 <- function(t, 
                                 t1,t2,
                                 c01,c12,c2)
{
  if(t<t1) res <- c01
  if(t>=t1 & t<=t2) res <- c12
  if(t>t2) res <- c2
  
  return(res)
}


SEIR.betacontrol.gamma <- function(t, x, parms)
{
  S <- x[1] 
  E <- x[2]
  I <- x[3]
  R <- x[4]
  C <- x[5] ## cumul number of cases
  D <- x[6] ## cumul number of disease-related deaths
  
  with(as.list(parms),
{ 
  dS <- mu - beta*control.transmission(t,tau,k)*S*I - mu*S 
  dE <- beta*control.transmission(t,tau,k)*S*I - sigma*E - mu*E
  dI <- sigma*E - (gamma+mu)*I
  dR <- (1-f)*gamma*I    # f = fatality rate
  
  dC <- sigma*E
  dD <- f*gamma*I
  
  out <- c(dS, dE, dI,dR,dC,dD) 
  list(out)
})
} 


SEIR.betacontrol2.gamma <- function(t, x, parms)
{
  S <- x[1] 
  E <- x[2]
  I <- x[3]
  R <- x[4]
  C <- x[5] ## cumul number of cases
  D <- x[6] ## cumul number of disease-related deaths
  
  with(as.list(parms),
       { 
         dS <- mu - beta*control.transmission2(t,t1,t2,c01,c12,c2)*S*I - mu*S 
         dE <- beta*control.transmission2(t,t1,t2,c01,c12,c2)*S*I - sigma*E - mu*E
         dI <- sigma*E - (gamma+mu)*I
         dR <- (1-f)*gamma*I    # f = fatality rate
         
         dC <- sigma*E
         dD <- f*gamma*I
         
         out <- c(dS, dE, dI,dR,dC,dD) 
         list(out)
       })
} 

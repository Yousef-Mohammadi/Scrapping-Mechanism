function dq = Cost_Fun( t,q,M,N,G,T )
global n_q
dq(1:2:2*n_q,1)=q(2:2:2*n_q);
dq(2:2:2*n_q,1)=M\(T+G-N*q(2:2:2*n_q));
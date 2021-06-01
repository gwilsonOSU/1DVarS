%
% perturbation test of tl_ws_brownLawler.m
%
clear

d50=200e-6;

verb=0;
for niter=1:100

% choose reasonable perturbations
frac_tl = 0.01;
myrand=@()2*(rand(1)-.5);
tl_d50=(d50-d50*(frac_tl*myrand()+1));

% compare TL perturbation to actual
ws=ws_brownLawler(d50);
ws_prime=ws_brownLawler(d50+tl_d50);
tl_ws(niter)=tl_ws_brownLawler(tl_d50,d50);

% calculate stats
tl_ws_true(niter) = ws_prime - ws;
percerr(niter)=max(abs((tl_ws_true(:,niter) - tl_ws(:,niter))./tl_ws_true(:,niter)))*100;
if(verb)
  disp(['true perturbation = ' num2str(ws_prime - ws)])
  disp(['TL   perturbation = ' num2str(tl_ws)])
  disp(['agreement to ' num2str(percerr(niter),3) ' %'])
end

end

% show all stats.  At first I was plotting percentage errors, but this
% turned out to be flawed as I was seeing a lot of division by small
% numbers.  Since then I changed to a scatterplot comparison, which is a
% much better indicator of validity
if(niter>1)
  clf
  plot(tl_ws_true,tl_ws,'.')
  hold on
  axis equal tight
  ax=axis;
  plot(ax([1 2]),ax([1 2]),'k--')
  hold off
  title(['ws comparison (should be 1-1)'])
  xlabel('true [m/s]')
  ylabel('predicted [m/s]')
end

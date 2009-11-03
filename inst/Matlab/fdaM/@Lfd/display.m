function display(Lfd)

%  Last modified 20 July 2006

nderiv = Lfd.nderiv;
fprintf(['NDERIV = ', num2str(nderiv),'\n']);
if nderiv > 0
    fprintf('\nWFD:');
    fprintf('\n\n-------------------');
    for ideriv=1:nderiv
        fprintf(['\n\nWFD(',num2str(ideriv-1),') fdPar object:\n'])
        display(Lfd.bwtcell{ideriv});
        fprintf('\n\n-------------------');
    end
end


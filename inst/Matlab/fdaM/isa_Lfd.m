function isaLfd = isa_Lfd(Lfdobj)
%  ISA_LFD  checks argument is either integer or functional data object.

%  last modified 2 December 2006

isaLfd = 0;

%  check class of LFDOBJ for Lfd

if strcmp(class(Lfdobj), 'lfd') || strcmp(class(Lfdobj), 'Lfd') 
    isaLfd = 1;
end

%  check class for being double, if so check that
%  it is a nonnegative integer

% if strcmp(class(Lfdobj), 'double') 
%     if floor(Lfdobj) == Lfdobj & Lfdobj >= 0
%         isaLfd = 1;
%     end
% end

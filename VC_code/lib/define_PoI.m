function PoI = define_PoI(MP,Nd,El)

PoI.cBM = median(find(strcmpi(Nd.name,'AM'))); % Center point of BM
PoI.zBM = find(strcmpi(Nd.name,'AP'));
PoI.zTM = find(strcmpi(Nd.name,'E1'));

PoI.pBM = median(find(strcmpi(Nd.name,'AM'))); % footage of DC
PoI.eHB = median(find(strcmp(El.name,'ANK'))); % HB
PoI.eOHC = median(find(strcmp(El.name,'OHC'))); % OHC
PoI.nRL = find(strcmp(Nd.name,'BB')); % RL
PoI.nTM = median(find(strcmp(Nd.name,'E1'))); % TM

PoI.aPC = find(strcmpi(Nd.name,'CC')); % apex of pillar cells
PoI.mPC = find(strcmpi(Nd.name,'C2')); % apex of pillar cells
PoI.rPC = find(strcmpi(Nd.name,'AP')); % apex of pillar cells
PoI.aDC = find(strcmpi(Nd.name,'DD')); % apex of Deiter's cells

zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
idx = knnsearch(MP.zz(:),zz(:));

fields = fieldnames(PoI);

for i = 1:length(fields)
    if length(PoI.(fields{i})) > 1
        PoI.(fields{i}) = PoI.(fields{i})(idx);
    end
end

end

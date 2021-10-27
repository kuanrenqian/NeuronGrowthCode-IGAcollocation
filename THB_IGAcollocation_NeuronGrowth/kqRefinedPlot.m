function [output] = kqRefinedPlot(input_cp,THBfinal,coll_X,coll_Y)

output = griddata(THBfinal(:,1),THBfinal(:,2),input_cp,coll_X,coll_Y);
imagesc(output);
axis square;
colorbar;
    
end
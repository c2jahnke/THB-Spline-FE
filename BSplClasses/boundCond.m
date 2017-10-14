classdef boundCond < handle
    properties
        west = 'Dirichlet';
        wVal = 0;
        eVal = 0;
        east = 'Neumann';
        mixed = true;
    end
    methods (Access = public)
    function obj = boundCond(west,east,eVal,wVal)
        %% (east,west,eVal,wVal)
        obj.east = east;
        obj.west = west;
        obj.eVal = eVal;
        obj.wVal = wVal;
        if(strcmp(east,west))
            obj.mixed = false;
        end
    end
    end

end
% convert a structure to a string
function str = struct2str(my_struct)
fields = sortrows(fieldnames(my_struct)); % avoid permutations
str = '';
for i=1:length(fields)
    fn = fields{i};
    fv = my_struct.(fn);
    switch class(fv)
        case 'char'
            str = sprintf('%s %s=%s',str,fn,fv);
        case 'double'
            if (round(fv)==fv)
                str = sprintf('%s %s=%d',str,fn,fv);
            else
                found = 0;
                for j=0:10
                    if (round(10^j*fv)==10^j*fv) && ~found
                        str = sprintf(sprintf('%%s %%s=%%.%df',j),str,fn,fv);
                        found = 1;
                    end
                end
            end
        otherwise
            error('unknown field class');
    end
end
str = str(2:end);

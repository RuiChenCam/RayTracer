function table = patterngen(funcHandle)
%This function generates the gain of an antenna from an input function
%handle and output a lookup table to be searched
%table in the form [theta phi gain]
    
    theta=-pi:0.01:pi;       
    phi=0:0.01:pi;
    table=zeros(length(theta)*length(phi),3);
    lenTheta=length(theta);
    lenPhi=length(phi);
    for a=1:lenTheta
        for b=1:lenPhi
            table((a-1)*lenPhi+b,:)=[theta(a),phi(b),funcHandle(theta(a),phi(b))];
        end
    end
end


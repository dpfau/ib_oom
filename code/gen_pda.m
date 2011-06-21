function [x z] = gen_pda(O,T,len,p)
% Generates a sequence from a probabilistic pushdown automata
% or should I call this a "Hidden Pushdown Model"...
% or "Hidden Markov Stack Model"

z = struct( 'state', {cell(len,1)}, 'stack', {cell(len,1)} );

nss = size(T,1)*size(T,2); % the number of states times 2 + the number of stack symbols (pop one, do nothing, or push that symbol)
x = zeros(len,1);
for i = 1:len
    z.state{i} = p.state;
    z.stack{i} = p.stack;
    if isempty( z.stack{i} )
        x(i) = sample( O( :, p.state, 1 ) );
        [a b] = ind2sub( sample( T( 1:nss ) ) );
        p.state = a;
        if b > 2
            p.stack = b - 2;
        end
    else
        x(i) = sample( O( :, p.state, p.stack(1) + 1 ) );
        [a b] = ind2sub( sample( T( p.stack(1)*nss + (1:nss) ) ) );
        p.state = a;
        if b == 2
            p.stack = p.stack(2:end);
        elseif b > 2
            p.stack = [ b - 2, p.stack ];
        end
    end
end
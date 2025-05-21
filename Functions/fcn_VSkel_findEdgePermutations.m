function edge_permutations = ...
    fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, varargin)

% Given the edges engaged by each vertex, and the total number of edges to
% consider, finds all the permutations of edge numbers that need to be
% checked. Any repeated permutations are removed.
%
% In 2D: the maximum number of permutations (rows) will be V*(NE-2) where V is the
% number of vertices, and NE is the number of edges
%
% In 3D: the maximum number of permutations (rows) will be 
%         V*(M choose 3)*(NE-M) 
% where V is the number of vertices, M is the number of edges in each
% vertex (can be 3 or more), and NE is the number of edges. Since M is
% often 3, and sometimes more, the above estimate will be V*(NE-3) or more.
%
% For example, consider a polytope with 20 vertices and 20 edges. Then vertex 1 of 20
% will "use" 2 edges to make that vertex - typically, these are in sequence
% (e.g. the vertex is between, for example, edge 4 and 5). This means for
% vertex 1, that 18 edges need to be tested against the 2 forming the
% vertex. 
% 
% Thus, edges 4 and 5 have to be tested against 1, 2, 3, 6, 7, ..., 20
%
% In this example, for all 2-edge combintations, there are 20*18 = 360
% permutations to test. This code calculates all permutations
%
% FORMAT:
%
% edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (fig_num))
%
% INPUTS:
%
%     cell_array_edges_in_vertices: a {V-by-1} cell array where each cell
%     contains a column vector of the edge IDs that each vertex contains.
%
%     NE: an integer representing the total number of edges.
%
%     (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     edge_permutations: a matrix of all permutations that need to be
%     tested, one in each row, containing no repeats and sorted,
%     column-wise, so that smalled edge IDs appear first.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_findEdgePermutations
%
% This function was written on 2025_05_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% REVISION HISTORY:
% 2025_05_18 by Sean Brennan
% -- first written by S. Brennan



% TO DO
% -- speed up by checking off which edges were checked previously, and not
% repeating the calculations again for these. For example, in a 4-sided 2D
% polytope, vertex 1 contains edges 41. If this hits edge 3, the result
% will be the same as vertex 4 which is edges 34 hitting edge 1. They all
% contain edges 134.

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
    flag_do_debug = 0; %     % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; %     % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS");
    MATLABFLAG_VSKEL_FLAG_DO_DEBUG = getenv("MATLABFLAG_VSKEL_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_VSKEL_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_VSKEL_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the cell_array_edges_in_vertices input
        assert(iscell(cell_array_edges_in_vertices))


        % Check the NE input
        fcn_DebugTools_checkInputsToFunctions(NE, 'positive_1column_of_integers',1);

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  3 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp; %#ok<NASGU>
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig = figure;
        fig_for_debug = fig.Number; %#ok<NASGU>
        flag_do_plot = 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the edges engaged by each vertex, and the total number of edges to
% consider, finds all the permutations of edge numbers that need to be
% checked. Any repeated permutations are removed.
%
% In 2D: the maximum number of permutations (rows) will be V*(NE-2) where V is the
% number of vertices, and NE is the number of edges
%
% In 3D: the maximum number of permutations (rows) will be 
%         V*(M choose 3)*(NE-M) 
% where V is the number of vertices, M is the number of edges in each
% vertex (can be 3 or more), and NE is the number of edges. Since M is
% often 3, and sometimes more, the above estimate will be V*(NE-3) or more.
%
% For example, consider a polytope with 20 vertices and 20 edges. Then vertex 1 of 20
% will "use" 2 edges to make that vertex - typically, these are in sequence
% (e.g. the vertex is between, for example, edge 4 and 5). This means for
% vertex 1, that 18 edges need to be tested against the 2 forming the
% vertex. 
% 
% Thus, edges 4 and 5 have to be tested against 1, 2, 3, 6, 7, ..., 20
%
% In this example, for all 2-edge combintations, there are 20*18 = 360
% permutations to test. This code calculates all permutations.
% NOTE: In 2D and 3D, this is MUCH smaller than N choose (D+1), e.g. N
% edges with 3 constraining edges. For example 20-choose-3 = 1140.
% (20*19*18)/(2*3)
% In the case of large vertex polytopes, say 400 vertices, this difference
% is substantial. For 400 vertices, the permutations to test is 159,200
% permutations to test. 400 choose 3 is 10,586,800, which is roughly 100
% times larger. The reason that the vertex count is so much smaller than
% the N-choose-M solution is that each vertex only allows certain other
% edges to be included, not just any combinations. 
%
% To solve, we permute along edges in vertices, which are numbered by
% vertex. 
% In 2d:
%      vertex 1 may contain edges 1 and 2, vertex 2 may contain 2 and 3,
%      vetext 3 may contain 7 and 8, etc. So the vertex listings would be:
%      [1 2; 2 3; 7 8; etc.]
%      so the allowable vertices would be:
%      [1 2 (3 to E); (1) 2 3; 2 3 (4 to E); (1 to 6) 7 8; 7 8 (9 to E);
%      etc.
%      Many of the above contain repeats, so these should be removed.
%      Observe that for vertices that have 2 edges, the number of
%      permutations per vertex is simply (E-2). These allowable
%      permutations can be generated by mapping the edge numbers (say, 7
%      and 8) to 1' and 2'. Thus, permuting through:
%      [1' 2' (3' to E')] generates all possible combinations in
%      prime-space. A substitution can then be used to find the original
%      valued permutations. For example, consider a polytope with 6 edges.
%      Then, if:
%
%      prime_permutations = [1 2 3; 1 2 4; 1 2 5; 1 2 6];
%
%      And the mapping from primes to normal edges is given by the
%      following, which is examining permutations engaging edges 2 and 3:
%
%      mapping_prime_to_normal = [2 3 1 4 5 6];
%
%      Then the "normal" permutations would be given by:
%
%      mapping_prime_to_normal(prime_permutations);
%
%     The above command says, basically, every place where there is a 1 in
%     prime_permutations, substitute the mapping number.
%      
%      This results in:
%         2     3     1
%         2     3     4
%         2     3     5
%         2     3     6
%
%     Note that the permutation template and mapping numbers can be
%     repeated to allow all permutations to be done at once. To do this,
%     the permutation_template needs to be offset by the number of
%     verticies:
%
%     prime_permutations_with_offsets = ...
%         [...
%         [1 2 3; 1 2 4; 1 2 5; 1 2 6];  % <--FIrst 2 permuted with 3:N
%         [1 2 3; 1 2 4; 1 2 5; 1 2 6] + 6;   % Offset by N
%         [1 2 3; 1 2 4; 1 2 5; 1 2 6] + 12]; % Offset by 2N        
%
%      mapping_prime_to_normal = [2 3 1 4 5 6, 4 5 1 2 3 6, 5 6 1 2 3 4];
%
%      mapping_prime_to_normal(prime_permutations_with_offsets)
%
%      Gives:
%
%         2     3     1
%         2     3     4
%         2     3     5
%         2     3     6
%         4     5     1
%         4     5     2
%         4     5     3
%         4     5     6
%         5     6     1
%         5     6     2
%         5     6     3
%         5     6     4
%
%
% In 3d:
%     Note that edges in 3D are called "faces", and their numbering is usually NOT
%     in sequence at a vertex, and a vertex may connect more than 3 faces.
%     Also, 4 edges are needed to define a point that is equidistant. Points
%     equidistant from 2 edges form a plane. Points equidistant from 3 form a
%     line. Each vertex connecting 3 faces thus defines 1 projection line.
%     with 4 faces, there are 4 choose 3 = 4 possible lines. With 5, there are
%     5 choose 3 or 10 possible lines.
% 
%     As well, for vertices with more than 3 faces, any combinations of face
%     tests should avoid the other faces also engaging that same vertex. So,
%     for example, imagine a vertex engaging faces 1, 5, 6, 9, and 12. If faces
%     1, 5, and 6 are chosen for testing, then 9 and 12 must be avoided in any
%     test combinations, since the equi-distant solution is the vertex itself.
% 
%     For vertices connecting M faces, the combinations of 3 faces have to be chosen
%     from this list, and then permutated while excluding others. In other
%     words, we would do (M choose 3) "generator" faces followed by (F - M)
%     other combinations with non-vertex faces. In the example above, this
%     results in:
%     [1 5 6 (excluding 9 and 12) - against all other edges]
%     [1 5 9 (excluding 6 and 12) - against all other edges]
%     etc.
%     Observe that one can map these non-sequential numbers to sequential ones:
%     1 -->  1'
%     5 -->  2'
%     6 -->  3'
%     9 -->  4'
%     12 --> 5'
%     all other numbers get mapped, in sequence, to 6', 7', 8', etc.
%     So, with the re-mapping, finding [1 5 6 (excluding 9 and 12) - against all other edges]
%     is the same as finding:
%     [1' 2' 3'] (skipping 4 and 5) tested against [6' to E']
%     We then "unmap" the results to get the permutations
% 
%     We can observe, that, in the 3D case, the number of permutations for a
%     vertex containing M edges will be (M choose 3)*(E-M)
% 
%     NOTE: Euler's formula for polyhedra: V-E+F=2, where V is the number of
%     vertices, E is the number of edges, and F is the number of faces
%
%     So, in the 3D case (with 5 faces, 4 at vertex 1, 3 at other vertices
%     (this is a square-bottomed pyramid, BTW):
% 
%                    V1
%                   /=\\
%                  /===\ \
%                 /=====\' \
%                /=======\'' \
%               /=========\ ' '\
%              /===========\''   \
%             /=============\ ' '  \
%            /===============\  F2'  \
%           /=======F1 =======\' ' ' ' \
%          /===================\' ' '  ' \
%         /=====================\' '   ' ' V3
%        /=======================\  '   ' /
%       /=========================\   ' /
%      /===========================\'  /
%     V5============================V2
%     
%     
%     Nfaces = 5;
%     prime_permutations_with_offsets_fourFaces = ...
%         [...
%         % Vertex 1 contains faces 1, 2, 3, 4, e.g. M = 4
%         [1 2 3 5];  % <--M choose 3 permuted with (M+1):N
%         [1 2 4 5];  % <--M choose 3 permuted with (M+1):N
%         [1 3 4 5];  % <--M choose 3 permuted with (M+1):N
%         [2 3 4 5];  % <--M choose 3 permuted with (M+1):N
%         ];
%     
%     prime_permutations_with_offsets_threeFaces = ...
%         [...
%         % Vertex 1 contains faces 1, 2, 3, 4, e.g. M = 4
%         [1 2 3 4];  % <--M choose 3 permuted with (M+1):N
%         [1 2 3 5];  % <--M choose 3 permuted with (M+1):N
%         ]; %#ok<NBRAK2>
%     
%     prime_permutations_with_offsets = [...
%         prime_permutations_with_offsets_fourFaces  + 0*Nfaces; ... % Vertex 1
%         prime_permutations_with_offsets_threeFaces + 1*Nfaces; ... % Vertex 2
%         prime_permutations_with_offsets_threeFaces + 2*Nfaces; ... % Vertex 3
%         prime_permutations_with_offsets_threeFaces + 3*Nfaces; ... % Vertex 4
%         prime_permutations_with_offsets_threeFaces + 4*Nfaces; ... % Vertex 5
%         ];
%     
%    mapping_prime_to_normal = [
%        1 2 3 4 5, ...  % Vertex 1
%        1 2 5 3 4, ...  % Vertex 2
%        2 3 5 1 4, ...  % Vertex 3
%        3 4 5 1 2, ...  % Vertex 4
%        4 1 5 2 3, ...  % Vertex 5
%    ];
%     
%    mapping_prime_to_normal(prime_permutations_with_offsets)
%    
%    % Gives:
%     1     2     3     5
%     1     2     4     5
%     1     3     4     5
%     2     3     4     5
%     1     2     5     3
%     1     2     5     4
%     2     3     5     1
%     2     3     5     4
%     3     4     5     1
%     3     4     5     2
%     4     1     5     2
%     4     1     5     3
%
% % As another example, if Nfaces = 8, and there are 4 faces at each
% % vertex:
% prime_permutations_with_offsets_fourFaces = ...
%     [...
%     % Vertex 1 contains faces 1, 2, 3, 4, e.g. M = 4
%     [1 2 3 5; 1 2 3 6; 1 2 3 7; 1 2 3 8];  % <--M choose 3 permuted with (M+1):N
%     [1 2 4 5; 1 2 4 6; 1 2 3 7; 1 2 4 8];  % <--M choose 3 permuted with (M+1):N
%     [1 3 4 5; 1 3 4 6; 1 3 4 7; 1 3 4 8];  % <--M choose 3 permuted with (M+1):N
%     [2 3 4 5; 2 3 4 6; 2 3 4 7; 2 3 4 8]];  % <--M choose 3 permuted with (M+1):N


% Is this 2D or 3D?
dimension_of_points = length(cell_array_edges_in_vertices{1});
if dimension_of_points<=1
    warning('on','backtrace');
    warning('Expecting a 2D or 3D point set. Dimension found: %.0f',dimension_of_points);
    error('Dimension 0 or 1 not coded yet.');
elseif dimension_of_points>=3
    dimension_of_points = 3;
end

Nvertices = length(cell_array_edges_in_vertices);

if 2==dimension_of_points
    %%%%
    % CREATE prime_permutations_with_offsets matrix:
    % Initialize prime_permutations array. Note that this is simply the
    % first array, shifted over/over/over
    prime_permutations_template = [1*ones(Nvertices-2,1) 2*ones(Nvertices-2,1) (3:Nvertices)'];
    prime_permutations_template_repeated = repmat(prime_permutations_template,Nvertices,1);

    % This next part is tricky. The goal is to create a matrix of offsets,
    % e.g. 0's for the first matrix, V's for the next matrix, 2V's for the
    % next matrix. To do this, we create an offset sequence that is the
    % same number of vertices.
    offset_sequence = (0:Nvertices-1)';

    % Now, we need this to have the same number of columns and rows as the
    % prime_permutations_template_repeated. This means the offset sequence
    % needs to be repeated (Nvertices - 2) times (e.g., the size of the
    % prime_permutations_template);
    Nrepeats = numel(prime_permutations_template);

    % Now make the offsets into columns where every row is 0, 1, 2, etc.
    % This requires a few transpose operations to work correctly
    offsets_repeated = repmat(offset_sequence,1,Nrepeats);
    offsets_shape_of_permulations = reshape(offsets_repeated',3,Nvertices*length(prime_permutations_template(:,1)))';
    prime_permutations_with_offsets = offsets_shape_of_permulations*Nvertices + prime_permutations_template_repeated;

    %%%%
    % CREATE mapping_prime_to_normal matrix:
    mapping_prime_to_normal_matrix = nan(Nvertices,Nvertices);
    edges_for_each_vertex = reshape(cell2mat(cell_array_edges_in_vertices),2,Nvertices)';
    
    % Fill in edges
    mapping_prime_to_normal_matrix(:,1:2) = edges_for_each_vertex;

    % Fill in remainders
    all_vertices = (1:Nvertices);
    for ith_vertex = 1:Nvertices
        mapping_prime_to_normal_matrix(ith_vertex,3:Nvertices) = setdiff(all_vertices,edges_for_each_vertex(ith_vertex,:));
    end
    mapping_prime_to_normal = reshape(mapping_prime_to_normal_matrix',1,Nvertices*Nvertices);

    edge_permutations_unsorted = mapping_prime_to_normal(prime_permutations_with_offsets);
    
elseif 3==dimension_of_points
    %%%%
    % CREATE prime_permutations_with_offsets matrix:
    
    % Find the number of faces in each vertex:
    NfacesEachVertex = zeros(Nvertices,1); % Initialize array
    for ith_vertex = 1:length(cell_array_edges_in_vertices)
        NfacesEachVertex(ith_vertex,1) = length(cell_array_edges_in_vertices{ith_vertex});
    end

    % Find which ones are unique
    templateSizesToPrecalculate = unique(NfacesEachVertex);

    URHERE
    Nvertices = length(cell_array_edges_in_vertices)

    % Initialize prime_permutations array. Note that this is simply the
    % first array, shifted over/over/over
    prime_permutations_template = [1*ones(Nvertices-2,1) 2*ones(Nvertices-2,1) (3:Nvertices)'];
    prime_permutations_template_repeated = repmat(prime_permutations_template,Nvertices,1);

    % This next part is tricky. The goal is to create a matrix of offsets,
    % e.g. 0's for the first matrix, V's for the next matrix, 2V's for the
    % next matrix. To do this, we create an offset sequence that is the
    % same number of vertices.
    offset_sequence = (0:Nvertices-1)';

    % Now, we need this to have the same number of columns and rows as the
    % prime_permutations_template_repeated. This means the offset sequence
    % needs to be repeated (Nvertices - 2) times (e.g., the size of the
    % prime_permutations_template);
    Nrepeats = numel(prime_permutations_template);

    % Now make the offsets into columns where every row is 0, 1, 2, etc.
    % This requires a few transpose operations to work correctly
    offsets_repeated = repmat(offset_sequence,1,Nrepeats);
    offsets_shape_of_permulations = reshape(offsets_repeated',3,Nvertices*length(prime_permutations_template(:,1)))';
    prime_permutations_with_offsets = offsets_shape_of_permulations*Nvertices + prime_permutations_template_repeated;

else
    warning('on','backtrace');
    warning('A vector was given that has dimension: %.0d, where 2D or 3D was expected',dimension_of_points);
    error('Function not yet coded for anything other than 2D or 3D');
end

% Sort by columns to make lowest indices first, then sort by entire rows
edge_permutations_sorted_columns = sort(edge_permutations_unsorted,2,"ascend");
edge_permutations_sorted = sortrows(edge_permutations_sorted_columns,"ascend");

[edge_permutations, IA] = unique(edge_permutations_sorted,'rows','legacy');

%% Plot results?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plot
    
    fprintf(1,'\nRESULTS of findEdgePermutations:\n');
    fprintf(1,'Dimension of search: %.0d\n',dimension_of_points);
    fprintf(1,'Number of edges (input NE): %.0d\n',NE);
    fprintf(1,'Unsorted permutations:\n')
    disp(edge_permutations_unsorted);
    fprintf(1,'Sorted permutations:\n')
    disp(edge_permutations_sorted);
    fprintf(1,'Tagged as unique:\n');
    for ith_permutation = 1:length(edge_permutations_sorted(:,1))
        for jth_permutation = 1:length(edge_permutations_sorted(1,:))
            fprintf('%.0f\t', edge_permutations_sorted(ith_permutation,jth_permutation))
        end
        if any(ith_permutation==IA)
            fprintf(1,' <--- \n');
        else
            fprintf(1,'\n');
        end
    end
    fprintf(1,'Final permutations (%.0d of them, from %.0d before removing repeats):\n', length(edge_permutations(:,1)), length(edge_permutations_sorted(:,1)));
    disp(edge_permutations);

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends INTERNAL_fcn_findUnitDirectionVectors



%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

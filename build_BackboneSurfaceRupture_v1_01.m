%% =========================================================================
%  BUILD_BACKBONESURFACERUPTURE
%
%  BACKBONE SURFACE RUPTURE (BSR)
%  Dynamic Backbone Graph + minimum-travel-time path selection
%
%  PURPOSE
%  -------
%  This script reconstructs a single ordered Backbone Surface Rupture (BSR)
%  from complex, multi-segment principal surface rupture traces.
%
%  The key difference from a feature-based Least-Cost Path is that the
%  original shapefile features are NOT treated as indivisible objects.
%
%  The observed principal-fault geometry is immutable. Before building the graph,
%  the code only removes consecutive duplicate vertices and zero-length
%  segments. It never merges features, moves tips, or replaces observed
%  polylines with simplified geometries.
%
%    1) preserves every cleaned original principal-fault feature;
%    2) identifies original tips and Douglas-Peucker significant vertices;
%    3) projects tips onto nearby features to define optional junctions;
%    4) splits graph edges at those positions while retaining the exact
%       original polyline geometry between consecutive graph nodes;
%    5) builds a topology above the immutable observed geometry;
%    6) requires each candidate to start and end on a TRACE edge;
%    7) forbids U-turn-like transitions greater than maximumAllowedTurn_deg;
%    8) computes minimum-travel-time paths between plausible endpoints;
%    9) selects the BSR lexicographically:
%
%         a) completion of the global rupture extent;
%         b) maximum trace coverage;
%         c) minimum number of artificial gaps;
%         d) minimum maximum gap;
%         e) minimum total gap;
%         f) minimum mean turn;
%         g) minimum backtracking;
%         h) minimum travel time.
%
%  This allows a path to enter a feature at an internal junction and use
%  only the required portion of that feature. It therefore avoids paths
%  that travel to a tip and then reverse direction along the same trace.
%
%  FINAL FIGURE STYLE
%  ------------------
%  Observed principal-fault traces: light grey, LineWidth 0.8.
%  Backbone Surface Rupture: black, LineWidth 1.5.
%  Artificial connections: red dashed, LineWidth 1.5.
%  Dynamic junctions remain in the graph/output but are not plotted.
%
%  INPUT
%  -----
%  Any polyline shapefile containing principal-fault traces.
%  Optionally select the principal-fault subset using any attribute field
%  and either a numeric or text value (e.g. Comp_rank = 1 or Class = 'Main').
%  If principalClassField is empty, all features are treated as principal.
%
%  OUTPUT
%  ------
%  BSR_RESULTS/IdE_<ID>/
%
%    BSR_IdE_<ID>.csv
%       ordered vertices, cumulative distance and x/L;
%
%    BSR_edges_IdE_<ID>.csv
%       ordered graph edges used by the selected BSR;
%
%    BSR_nodes_IdE_<ID>.csv
%       nodes of the dynamic graph;
%
%    BSR_all_edges_IdE_<ID>.csv
%       all trace and artificial-gap edges;
%
%    BSR_candidates_IdE_<ID>.csv
%       all valid candidate paths and selection diagnostics;
%
%    BSR_IdE_<ID>.png / .fig
%       diagnostic figure;
%
%    BSR_selection_IdE_<ID>.txt
%       explanation of the hierarchical selection;
%
%  Global summary:
%    BSR_DYNAMIC/BSR_summary.csv
%
%  REQUIREMENTS
%  ------------
%  Mapping Toolbox: shaperead.
%
%  Coordinate note:
%  The script uses a local equirectangular projection centred on each
%  event. For rupture extents of a few tens of kilometres, this is normally
%  sufficient for the present geometric analysis. It can be replaced by a
%  UTM conversion if desired.
% =========================================================================

clear;
clc;
close all;

%% =========================================================================
% USER SETTINGS
% =========================================================================

ruptureShp = fullfile( ...
    'SURE-main', ...
    'SURE2.0_ruptures', ...
    'SURE2.0_ruptures.shp');

% Field containing the event identifier. Event identifiers are expected to
% be numeric in this MATLAB reference implementation (e.g. SURE field 'IdE').
eventField = 'IdE';

% OPTIONAL PRINCIPAL-FAULT FILTER.
% Leave principalClassField empty ('') when the input shapefile already
% contains only the principal-fault traces. Otherwise specify any attribute
% field and the value identifying the principal traces.
%
% Examples:
%   SURE 2.0:          principalClassField = 'Comp_rank';
%                      principalClassValue = 1;
%   Generic dataset:   principalClassField = 'Class';
%                      principalClassValue = 'Main';
%   Pre-filtered layer: principalClassField = '';
%                      principalClassValue = [];
principalClassField = 'Comp_rank';
principalClassValue = 1;

% Empty means all events represented by the selected principal-fault traces.
% Alternatively provide a numeric vector, e.g. [19501214 19590818].
eventIDs = [];

% These lists are OPTIONAL and are used only to add Mw and kinematics to the
% global summary. They no longer control which events are processed.
normalListFile     = 'list_Normal.txt';
reverseListFile    = 'list_Reverse.txt';
strikeslipListFile = 'list_StrikeSlip.txt';

% -------------------------------------------------------------------------
% DYNAMIC JUNCTIONS AND FEATURE SPLITTING
% -------------------------------------------------------------------------

% A tip closer than this value to the interior of another feature creates
% a junction and splits the receiving feature.
% Increase this if visually plausible tip-to-feature links are missed.
junctionTolerance_m = 300;

% A projected junction must be at least this far from both ends of the
% receiving feature. This avoids unnecessary splits almost coincident with
% existing tips.
minimumSplitDistance_m = 25;

% Split positions closer than this value are merged.
splitMergeTolerance_m = 5;

% -------------------------------------------------------------------------
% IMMUTABLE OBSERVED GEOMETRY
% -------------------------------------------------------------------------

% Consecutive vertices closer than this value are removed. This is the only
% operation permitted to change the stored vertex sequence.
duplicateVertexTolerance_m = 0.05;

% Features shorter than this value after removal of duplicate vertices are
% discarded as degenerate.
minimumCleanFeatureLength_m = 1;

% Graph endpoints are clustered only when they are effectively coincident.
% This tolerance must remain small: it is not a geometric tip snap.
nodeSnapTolerance_m = 1.0;

% The following legacy preprocessing operations are deliberately disabled
% in V2. They remain absent from the processing call below.
preprocessingTipSnapTolerance_m = 0;
mergeFeatureTipTolerance_m = 0;
mergeFeatureMaximumTurn_deg = 0;
mergeMaximumTortuosity = Inf;
rejectSelfIntersectingMerges = false;
maximumMergePasses = 0;

% -------------------------------------------------------------------------
% DOUGLAS-PEUCKER SIGNIFICANT EXIT NODES
% -------------------------------------------------------------------------

% Internal vertices retained by Douglas-Peucker become graph nodes from
% which the path may leave a feature before its natural tip. The original
% feature geometry is still retained on every trace edge and in the output.
minimumSimplificationTolerance_m = 50;
maximumSimplificationTolerance_m = 500;
simplificationLengthFraction = 0.002;

% To keep the connection pool small, each significant exit node is linked
% only to the nearest original tips within maxGap_m.
maximumConnectionsPerExit = 50;

% -------------------------------------------------------------------------
% ARTIFICIAL GAPS
% -------------------------------------------------------------------------

% Maximum straight-line gap allowed between graph nodes.
maxGap_m = 15000;

% Beyond this gap length, the graph cost increases more rapidly.
preferredGap_m = 10000;

% Gap connections can be:
%   'terminalToAll'  - at least one endpoint must be a terminal trace node;
%   'terminalOnly'   - both endpoints must be terminal trace nodes;
%   'all'            - every node pair can be linked (can be expensive).
%   'significantToTerminal' - Douglas-Peucker node to original feature tip.
gapConnectionMode = 'significantToTerminal';

% Do not create artificial gaps shorter than this value because such nodes
% should normally already have been snapped together.
minimumArtificialGap_m = 2*nodeSnapTolerance_m;

% Final geometric continuity check. If two consecutive selected
% edges do not meet within this tolerance, insert an explicit GAP
% before recomputing metrics and assembling the exported BSR.
assemblyContinuityTolerance_m = 1.0;

% -------------------------------------------------------------------------
% TURN CONSTRAINT
% -------------------------------------------------------------------------

% A transition with a direction change larger than this value is excluded
% from the state graph. This prevents U-turns and strong reversals.
maximumAllowedTurn_deg = 140;

% -------------------------------------------------------------------------
% ENDPOINT CANDIDATES
% -------------------------------------------------------------------------

% Endpoint search is restricted to opposite sides of the rupture.
% At each level, the code selects the N original tips with the lowest
% progress and the N tips with the highest progress along the global PCA
% axis, then evaluates only start-side to end-side combinations. This
% avoids testing all possible terminal pairs, including internal tips.
%
% The search expands only if no candidate reaches the required global
% extent. With levels [10 20 30 50], at most 100, 400, 900 and 2500 pairs
% are tested at each successive level (before duplicate removal).
extremeTipsPerSideLevels = [10 20 30 50];
maxCandidatePairs = Inf;

% Only original feature tips are used as global start/end candidates.
% Dynamic interior junctions are graph nodes but not global endpoints.

% -------------------------------------------------------------------------
% GLOBAL-ENDPOINT COMPLETION
% -------------------------------------------------------------------------

% Candidate endpoint span is measured along the principal axis of all
% observed rank-1 geometry. It is considered before trace coverage, so a
% path cannot stop before the last principal rank-1 piece merely to save a
% connection.
globalExtentTolerance = 0.005;
minimumGlobalExtentForCompletion = 0.98;

% -------------------------------------------------------------------------
% LEXICOGRAPHIC SELECTION
% -------------------------------------------------------------------------

% Candidates within this amount of the maximum coverage pass to the next
% criterion. Example: max = 1.0 and tolerance = 0.005 -> keep >= 0.995.
coverageTolerance = 0.005;

nGapTolerance       = 0;
maxGapTolerance_m   = 1;
totalGapTolerance_m = 1;
meanTurnTolerance   = 0.1;
backtrackTolerance  = 1e-4;
pathCostTolerance   = 1e-6;

% -------------------------------------------------------------------------
% LENGTH AND X/L
% -------------------------------------------------------------------------

% Artificial gaps always contribute to the total BSR length and x/L.
% A trace-only cumulative distance is retained as a diagnostic output.

% -------------------------------------------------------------------------
% CANDIDATE FIGURES
% -------------------------------------------------------------------------

numberCandidateFigures = 3;  % Inf saves all candidates

% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------

outDir = 'BSR_RESULTS_FINAL';

% -------------------------------------------------------------------------
% MINIMUM-TRAVEL-TIME SCENARIOS
% -------------------------------------------------------------------------
%
% Travel speed on observed rank-1 TRACE edges is the reference speed:
%
%     v_trace = 1
%
% JUNCTION and GAP edges are travelled more slowly.  A multiplier K means:
%
%     v_off_trace = 1/K
%     travel time outside TRACE = K * geometric length
%
% Three scenarios are evaluated automatically.  Turn is used only as a
% hard admissibility constraint; it is not included in travel time.
%
offTraceTimeMultipliers = 5;

travelTimeScenarios = repmat(struct( ...
    'traceSpeed',1.0, ...
    'offTraceSpeed',NaN, ...
    'offTraceMultiplier',NaN),numel(offTraceTimeMultipliers),1);

for iScenario = 1:numel(offTraceTimeMultipliers)
    travelTimeScenarios(iScenario).traceSpeed = 1.0;
    travelTimeScenarios(iScenario).offTraceMultiplier = ...
        offTraceTimeMultipliers(iScenario);
    travelTimeScenarios(iScenario).offTraceSpeed = ...
        1/offTraceTimeMultipliers(iScenario);
end

% Save one comparison figure per multiplier for every event.
saveTravelTimeScenarioFigures = true;

%% =========================================================================
% PREPARATION
% =========================================================================

if ~exist(outDir,'dir')
    mkdir(outDir);
end

fprintf('Reading %s ...\n',ruptureShp);
Rall = shaperead(ruptureShp);

% Event identifiers.
if ~isfield(Rall,eventField)
    error('Event field "%s" not found in the input shapefile.',eventField);
end
idAll = getNumericField(Rall,eventField);

% Principal-fault selection. The field is optional and the requested value
% can be numeric or text. Matching is case-insensitive for text and uses a
% small numerical tolerance for numeric values, mirroring the QGIS tool.
if strlength(strtrim(string(principalClassField)))==0
    principalMask = true(numel(Rall),1);
    fprintf('Principal-fault filter: none | using all input features.\n');
else
    if ~isfield(Rall,principalClassField)
        error('Principal-fault classification field "%s" not found.', ...
            principalClassField);
    end
    principalMask = matchAttributeValue( ...
        Rall,principalClassField,principalClassValue);
    fprintf('Principal-fault filter: %s = %s | %d input features selected.\n', ...
        principalClassField,char(string(principalClassValue)),sum(principalMask));
end

if ~any(principalMask)
    error('No input features match the selected principal-fault filter.');
end

eventInfo = readEventLists(normalListFile,reverseListFile,strikeslipListFile);

if isempty(eventIDs)
    eventIDs = unique(idAll(principalMask & isfinite(idAll)));
else
    eventIDs = eventIDs(:);
end

fprintf('Events to process: %d\n',numel(eventIDs));

summaryRows = {};

%% =========================================================================
% EVENT LOOP
% =========================================================================

for ie = 1:numel(eventIDs)

    eventID = eventIDs(ie);

    fprintf('\n[%d/%d] Event %d\n',ie,numel(eventIDs),eventID);

    idx = idAll == eventID & principalMask;
    PFraw = Rall(idx);

    if isempty(PFraw)
        warning('Event %d: no principal-fault features.',eventID);
        continue;
    end

    %% --------------------------------------------------------------------
    % Extract original polyline parts and convert to local metric coordinates.
    % ---------------------------------------------------------------------

    original = extractPolylineParts(PFraw);

    if isempty(original)
        warning('Event %d: no valid geometry.',eventID);
        continue;
    end

    allLon = vertcat(original.lon);
    allLat = vertcat(original.lat);

    lon0 = mean(allLon,'omitnan');
    lat0 = mean(allLat,'omitnan');

    for i = 1:numel(original)

        [original(i).x,original(i).y] = lonlat2local( ...
            original(i).lon,original(i).lat,lon0,lat0);

        original(i).cum_m = polylineCumulativeDistance( ...
            original(i).x,original(i).y);

        original(i).length_m = original(i).cum_m(end);
        original(i).startXY = [original(i).x(1),original(i).y(1)];
        original(i).endXY   = [original(i).x(end),original(i).y(end)];
    end

    valid = [original.length_m] > 0;
    original = original(valid);

    if isempty(original)
        warning('Event %d: all principal-fault parts are degenerate.',eventID);
        continue;
    end

    nPartsBeforePreprocessing = numel(original);

    [original,preprocessingStats] = cleanObservedRank1Geometry( ...
        original, ...
        duplicateVertexTolerance_m, ...
        minimumCleanFeatureLength_m);

    if isempty(original)
        warning('Event %d: preprocessing removed all principal-fault parts.',eventID);
        continue;
    end

    totalObservedLength_m = sum([original.length_m]);

    fprintf('  Original principal-fault parts: %d\n',nPartsBeforePreprocessing);
    fprintf('  Cleaned immutable principal-fault parts: %d\n',numel(original));
    fprintf('  Removed duplicate vertices: %d\n', ...
        preprocessingStats.RemovedVertices);
    fprintf('  Geometrically moved tips: 0\n');
    fprintf('  Merged feature pairs: 0\n');
    fprintf('  Total observed principal-fault length: %.2f km\n', ...
        totalObservedLength_m/1000);

    %% --------------------------------------------------------------------
    % Build dynamic split positions from tip-to-feature projections.
    % ---------------------------------------------------------------------

    [splitPositions,junctionRecords] = findDynamicJunctions( ...
        original, ...
        junctionTolerance_m, ...
        minimumSplitDistance_m, ...
        splitMergeTolerance_m);

    fprintf('  Dynamic interior junctions: %d\n', ...
        height(junctionRecords));

    % Add Douglas-Peucker vertices as internal split positions. These are
    % significant graph nodes and possible early-exit locations.
    [splitPositions,significantVertexRecords] = ...
        addDouglasPeuckerSplitPositions( ...
            original,splitPositions,splitMergeTolerance_m, ...
            minimumSimplificationTolerance_m, ...
            maximumSimplificationTolerance_m, ...
            simplificationLengthFraction);

    fprintf('  Douglas-Peucker significant vertices: %d\n', ...
        height(significantVertexRecords));

    %% --------------------------------------------------------------------
    % Split each original feature at all dynamic junction positions.
    % ---------------------------------------------------------------------

    traceEdgesRaw = splitOriginalFeatures( ...
        original,splitPositions,splitMergeTolerance_m);

    fprintf('  Trace subsegments after splitting: %d\n', ...
        numel(traceEdgesRaw));

    %% --------------------------------------------------------------------
    % Build and snap graph nodes.
    % ---------------------------------------------------------------------

    [nodes,traceEdges] = buildSnappedTraceGraph( ...
        traceEdgesRaw,nodeSnapTolerance_m,original);

    nodes.IsSignificantExit = identifySignificantExitNodes( ...
        nodes,significantVertexRecords,2*nodeSnapTolerance_m);

    % Original tips are also valid exit nodes.
    nodes.IsSignificantExit = ...
        nodes.IsSignificantExit | nodes.IsOriginalTip;

    traceDegree = accumarray( ...
        [traceEdges.Node1;traceEdges.Node2], ...
        1, ...
        [height(nodes),1]);

    nodes.TraceDegree = traceDegree;
    nodes.IsTerminal = traceDegree == 1;
    nodes.IsJunction = traceDegree >= 3;

    fprintf('  Graph nodes: %d\n',height(nodes));
    fprintf('  Terminal nodes: %d\n',sum(nodes.IsTerminal));
    fprintf('  Branch/junction nodes: %d\n',sum(nodes.IsJunction));
    fprintf('  Significant exit nodes: %d\n',sum(nodes.IsSignificantExit));

    %% --------------------------------------------------------------------
    % Add explicit tip-to-feature junction connectors.
    %
    % A projected tip and the new internal split node are not necessarily
    % coincident. The short connector below makes the dynamic junction
    % topologically usable by the graph without forcing the path to travel
    % to the original end of the receiving feature.
    % ---------------------------------------------------------------------

    junctionEdges = buildJunctionConnectorEdges( ...
        nodes,junctionRecords,nodeSnapTolerance_m);

    fprintf('  Junction connector edges: %d\n',height(junctionEdges));

    %% --------------------------------------------------------------------
    % Add artificial gap edges.
    % ---------------------------------------------------------------------

    existingEdges = concatenateEdgeTables(traceEdges,junctionEdges);

    % Recompute graph degree after adding the explicit junction connectors.
    graphDegree = accumarray( ...
        [existingEdges.Node1;existingEdges.Node2], ...
        1, ...
        [height(nodes),1]);

    nodes.GraphDegree = graphDegree;
    nodes.IsTerminal = graphDegree == 1;
    nodes.IsJunction = graphDegree >= 3;

    gapEdges = buildArtificialGapEdges( ...
        nodes, ...
        existingEdges, ...
        maxGap_m, ...
        minimumArtificialGap_m, ...
        gapConnectionMode, ...
        maximumConnectionsPerExit);

    allEdges = concatenateEdgeTables(existingEdges,gapEdges);

    fprintf('  Artificial gap edges available: %d\n',height(gapEdges));

    %% --------------------------------------------------------------------
    % Global endpoint span along the principal rupture axis.
    % ---------------------------------------------------------------------

    [globalOrigin,globalAxis,globalNodeProgress,globalObservedExtent_m] = ...
        computeGlobalEndpointProgress(original,nodes);

    %% --------------------------------------------------------------------
    % Determine candidate start/end nodes from original feature tips.
    %
    % Only pairs connecting opposite ends of the rupture are evaluated.
    % The N lowest-progress tips form the start-side pool and the N
    % highest-progress tips form the end-side pool. N is increased only
    % when no path reaches the required global extent.
    % ---------------------------------------------------------------------

    tipNodes = identifyOriginalTipNodes( ...
        original,nodes,nodeSnapTolerance_m);

    candidateTemplate = emptyCandidate();
    candidates = repmat(candidateTemplate,0,1);
    nCand = 0;

    searchedPairs = zeros(0,2);
    globalPairCounter = 0;
    completePathFound = false;

    for iLevel = 1:numel(extremeTipsPerSideLevels)

        nExtremeTipsPerSide = extremeTipsPerSideLevels(iLevel);

        [candidateTips,candidatePairs,startSideTips,endSideTips] = ...
            chooseOppositeExtremeTipPairs( ...
                tipNodes, ...
                globalNodeProgress, ...
                nExtremeTipsPerSide, ...
                maxCandidatePairs);

        % Remove endpoint pairs already evaluated at a previous level.
        pairMatrix = sort([candidatePairs.StartNode candidatePairs.EndNode],2);
        if ~isempty(searchedPairs)
            isNewPair = ~ismember(pairMatrix,searchedPairs,'rows');
            candidatePairs = candidatePairs(isNewPair,:);
            pairMatrix = pairMatrix(isNewPair,:);
        end

        fprintf(['  Endpoint-search level %d per side: ', ...
                 '%d start-side tips, %d end-side tips, %d new pairs\n'], ...
            nExtremeTipsPerSide, ...
            height(startSideTips), ...
            height(endSideTips), ...
            height(candidatePairs));

        fprintf('    candidate endpoint nodes represented: %d\n', ...
            height(candidateTips));

        if isempty(candidatePairs)
            continue;
        end

        levelCandidateStart = nCand + 1;

        for ipair = 1:height(candidatePairs)

            startNode = candidatePairs.StartNode(ipair);
            endNode   = candidatePairs.EndNode(ipair);
            globalPairCounter = globalPairCounter + 1;

            for iw = 1:numel(travelTimeScenarios)

                candidate = solveDynamicStatePath( ...
                    nodes, ...
                    allEdges, ...
                    startNode, ...
                    endNode, ...
                    travelTimeScenarios(iw), ...
                    preferredGap_m, ...
                    maximumAllowedTurn_deg);

                if ~candidate.isValid
                    continue;
                end

                candidate.PairIndex = globalPairCounter;
                candidate.WeightIndex = iw;
                candidate.StartNode = startNode;
                candidate.EndNode = endNode;

                candidate = evaluateDynamicCandidate( ...
                    candidate, ...
                    nodes, ...
                    allEdges, ...
                    totalObservedLength_m);

                candidate.GlobalExtentCoverage = min(1,abs( ...
                    globalNodeProgress(endNode)- ...
                    globalNodeProgress(startNode))/ ...
                    max(globalObservedExtent_m,eps));

                candidate.ReachesGlobalEndpoints = ...
                    candidate.GlobalExtentCoverage >= ...
                    minimumGlobalExtentForCompletion;

                candidate = orderfields(candidate,candidateTemplate);

                nCand = nCand+1;
                candidates(nCand,1) = candidate;
            end
        end

        searchedPairs = unique([searchedPairs;pairMatrix],'rows','stable');

        if nCand >= levelCandidateStart
            levelCandidates = candidates(levelCandidateStart:nCand);
            maxExtentLevel = max([levelCandidates.GlobalExtentCoverage]);
            completePathFound = any( ...
                [levelCandidates.ReachesGlobalEndpoints]);

            fprintf('    valid paths: %d | maximum extent: %.2f %%', ...
                numel(levelCandidates),100*maxExtentLevel);

            if completePathFound
                fprintf(' | complete path found\n');
            else
                fprintf(' | expanding extreme-tip search\n');
            end
        else
            fprintf('    valid paths: 0 | expanding extreme-tip search\n');
        end

        if completePathFound
            break;
        end
    end

    fprintf('  Endpoint pairs evaluated in total: %d\n', ...
        size(searchedPairs,1));

    if isempty(candidates)
        warning('Event %d: no valid dynamic path found.',eventID);
        continue;
    end

    candidates = removeDuplicateDynamicCandidates(candidates);

    selectionSettings = struct( ...
        'minimumGlobalExtentForCompletion', ...
            minimumGlobalExtentForCompletion, ...
        'globalExtentTolerance',globalExtentTolerance, ...
        'coverageTolerance',coverageTolerance, ...
        'nGapTolerance',nGapTolerance, ...
        'maxGapTolerance_m',maxGapTolerance_m, ...
        'totalGapTolerance_m',totalGapTolerance_m, ...
        'meanTurnTolerance',meanTurnTolerance, ...
        'backtrackTolerance',backtrackTolerance, ...
        'pathCostTolerance',pathCostTolerance);

    [selectedIndex,selectionLog,candidates] = ...
        selectCandidateLexicographically( ...
            candidates,selectionSettings);

    best = candidates(selectedIndex);

    % --------------------------------------------------------------------
    % Normalize the selected geometry for final output.
    % JUNCTION edges are internal graph objects. In the exported BSR they
    % must be represented as GAP edges whenever they span a finite distance.
    % Therefore the final BSR contains only TRACE and GAP geometries.
    % --------------------------------------------------------------------
    [allEdges,best.DirectedEdgeStates,nConvertedJunctions] = ...
        convertSelectedJunctionsToGaps( ...
            allEdges,best.DirectedEdgeStates, ...
            assemblyContinuityTolerance_m);

    % Final geometric continuity check. This does not alter path selection.
    % It inserts an explicit GAP wherever consecutive selected geometries
    % remain spatially discontinuous.
    [allEdges,best.DirectedEdgeStates,nInsertedAssemblyGaps] = ...
        insertMissingAssemblyGaps( ...
            nodes,allEdges,best.DirectedEdgeStates, ...
            assemblyContinuityTolerance_m);

    if nConvertedJunctions > 0 || nInsertedAssemblyGaps > 0
        best = evaluateDynamicCandidate( ...
            best,nodes,allEdges,totalObservedLength_m);
        candidates(selectedIndex) = best;

        if nConvertedJunctions > 0
            fprintf(['  Final output converted %d dynamic junction(s) ', ...
                'to artificial connection(s).\n'],nConvertedJunctions);
        end

        if nInsertedAssemblyGaps > 0
            fprintf('  Final continuity check inserted %d artificial connection(s).\n', ...
                nInsertedAssemblyGaps);
        end
    end

    fprintf('  Selected BSR:\n');
    fprintf('    global extent: %.2f %%\n',100*best.GlobalExtentCoverage);
    fprintf('    reaches global endpoints: %d\n',best.ReachesGlobalEndpoints);
    fprintf('    trace coverage: %.2f %%\n',100*best.Coverage);
    fprintf('    trace edges: %d\n',best.NumberTraceEdges);
    fprintf('    artificial gaps: %d\n',best.NumberArtificialGaps);
    fprintf('    maximum gap: %.1f m\n',best.MaximumGap_m);
    fprintf('    total gap: %.1f m\n',best.TotalGap_m);
    fprintf('    mean turn: %.1f deg\n',best.MeanTurn_deg);
    fprintf('    backtracking: %.5f\n',best.BacktrackingFraction);

    %% --------------------------------------------------------------------
    % Assemble the selected ordered BSR geometry.
    % ---------------------------------------------------------------------

    BSR = assembleDynamicBSR( ...
        nodes,allEdges,best.DirectedEdgeStates,lon0,lat0);

    fractionOnRank1 = BSR.TraceLength_m/max(BSR.LengthWithGaps_m,eps);
    [reconstructionConfidence,geometryFlag] = ...
        assignReconstructionConfidence( ...
            best.GlobalExtentCoverage, ...
            best.Coverage, ...
            fractionOnRank1, ...
            best.MaximumGap_m, ...
            best.BacktrackingFraction);

    best.ReconstructionConfidence = reconstructionConfidence;
    best.GeometryFlag = geometryFlag;
    best.FractionOnRank1 = fractionOnRank1;

    eventDir = fullfile(outDir,sprintf('IdE_%d',eventID));

    if ~exist(eventDir,'dir')
        mkdir(eventDir);
    end

    %% --------------------------------------------------------------------
    % Write BSR vertex table.
    % ---------------------------------------------------------------------

    nV = numel(BSR.X);

    BSRTable = table( ...
        repmat(eventID,nV,1), ...
        (1:nV)', ...
        BSR.X(:),BSR.Y(:), ...
        BSR.Lon(:),BSR.Lat(:), ...
        BSR.CumDistanceUsed_m(:), ...
        BSR.XoverL(:), ...
        BSR.CumDistanceTraceOnly_m(:), ...
        BSR.XoverL_TraceOnly(:), ...
        BSR.EdgeID(:), ...
        BSR.EdgeType(:), ...
        BSR.OriginalPartIndex(:), ...
        BSR.IdS(:), ...
        'VariableNames',{ ...
        'IdE','VertexID','X_m','Y_m','Longitude','Latitude', ...
        'AlongBSRDistance_m','XoverL', ...
        'AlongTraceDistance_m','XoverL_TraceOnly', ...
        'EdgeID','EdgeType','OriginalPartIndex','IdS'});

    writetable(BSRTable,fullfile( ...
        eventDir,sprintf('BSR_IdE_%d.csv',eventID)));

    %% --------------------------------------------------------------------
    % Write ordered selected-edge table.
    % ---------------------------------------------------------------------

    SelectedEdgeTable = buildSelectedEdgeTable( ...
        eventID,best,allEdges);

    writetable(SelectedEdgeTable,fullfile( ...
        eventDir,sprintf('BSR_edges_IdE_%d.csv',eventID)));

    %% --------------------------------------------------------------------
    % Write graph node and all-edge tables.
    % ---------------------------------------------------------------------

    nodesOut = nodes;
    nodesOut.IdE = repmat(eventID,height(nodesOut),1);
    nodesOut = movevars(nodesOut,'IdE','Before',1);

    writetable(nodesOut,fullfile( ...
        eventDir,sprintf('BSR_nodes_IdE_%d.csv',eventID)));

    allEdgesOut = removevars(allEdges,{'XGeom','YGeom'});
    allEdgesOut.IdE = repmat(eventID,height(allEdgesOut),1);
    allEdgesOut = movevars(allEdgesOut,'IdE','Before',1);

    writetable(allEdgesOut,fullfile( ...
        eventDir,sprintf('BSR_all_edges_IdE_%d.csv',eventID)));

    %% --------------------------------------------------------------------
    % Write candidate table and selection explanation.
    % ---------------------------------------------------------------------

    CandidateTable = candidatesToDynamicTable(candidates,eventID);

    writetable(CandidateTable,fullfile( ...
        eventDir,sprintf('BSR_candidates_IdE_%d.csv',eventID)));

    writeSelectionExplanation( ...
        fullfile(eventDir,sprintf('BSR_selection_IdE_%d.txt',eventID)), ...
        eventID,best,selectionLog);

    %% --------------------------------------------------------------------
    % Figures.
    % ---------------------------------------------------------------------

    makeBSRDiagnosticFigure( ...
        original,nodes,allEdges,best,BSR,eventID,eventDir);

    if saveTravelTimeScenarioFigures
        makeTravelTimeScenarioFigures( ...
            original,nodes,allEdges,candidates, ...
            travelTimeScenarios,eventID,eventDir, ...
            lon0,lat0,selectionSettings);
    end

    makeDynamicCandidateFigures( ...
        original,nodes,allEdges,candidates,eventID,eventDir, ...
        lon0,lat0,numberCandidateFigures);

    %% --------------------------------------------------------------------
    % Save MAT.
    % ---------------------------------------------------------------------

    save(fullfile(eventDir,sprintf('BSR_IdE_%d.mat',eventID)), ...
        'BSR','best','candidates','selectionLog', ...
        'original','nodes','traceEdges','junctionEdges','gapEdges','allEdges', ...
        'junctionRecords','significantVertexRecords','candidateTips','candidatePairs', ...
        'lon0','lat0');

    %% --------------------------------------------------------------------
    % Global summary row.
    % ---------------------------------------------------------------------

    infoRow = eventInfo(eventInfo.IdE == eventID,:);

    mw = NaN;
    kin = "";

    if ~isempty(infoRow)
        mw = infoRow.Mw(1);
        kin = infoRow.Kinematics(1);
    end

    summaryRows(end+1,:) = { ...
        eventID,mw,kin, ...
        numel(original),height(traceEdges),height(nodes), ...
        height(junctionRecords),totalObservedLength_m, ...
        best.Coverage,best.NumberTraceEdges,best.NumberJunctionEdges, ...
        best.NumberArtificialGaps,best.MaximumGap_m, ...
        best.MeanGap_m,best.TotalGap_m,best.MeanTurn_deg, ...
        best.BacktrackingFraction,best.StartNode,best.EndNode, ...
        best.WeightIndex,best.CandidateID,best.CandidateRank, ...
        BSR.TraceLength_m,BSR.LengthWithGaps_m, ...
        BSR.TotalLengthUsed_m,best.GlobalExtentCoverage, ...
        best.ReachesGlobalEndpoints,fractionOnRank1, ...
        reconstructionConfidence,geometryFlag};
end

%% =========================================================================
% GLOBAL SUMMARY
% =========================================================================

if ~isempty(summaryRows)

    Summary = cell2table(summaryRows,'VariableNames',{ ...
        'IdE','Mw','Kinematics', ...
        'NOriginalRank1Parts','NTraceSubsegments','NGraphNodes', ...
        'NDynamicJunctions','TotalObservedRank1Length_m', ...
        'Coverage','NTraceEdgesInBSR','NJunctionEdgesInBSR','NumberArtificialGaps', ...
        'MaximumGap_m','MeanGap_m','TotalGap_m','MeanTurn_deg', ...
        'BacktrackingFraction','StartNode','EndNode','TravelTimeScenario', ...
        'CandidateID','CandidateRank','BSR_TraceLength_m', ...
        'BSR_LengthWithGaps_m','BSR_LengthUsed_m', ...
        'GlobalExtentCoverage','ReachesGlobalEndpoints', ...
        'FractionOnRank1','ReconstructionConfidence','GeometryFlag'});

    writetable(Summary,fullfile(outDir,'BSR_summary.csv'));
end

fprintf('\nBSR processing completed.\n');

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function eventInfo = readEventLists(normalFile,reverseFile,strikeslipFile)

    eventInfo = table();

    files = {normalFile, reverseFile, strikeslipFile};
    labels = ["Normal", "Reverse", "StrikeSlip"];

    for iFile = 1:numel(files)
        filename = files{iFile};
        if exist(filename,'file')
            A = readmatrix(filename);
            if ~isempty(A)
                T = table(A(:,1),A(:,2),repmat(labels(iFile),size(A,1),1), ...
                    'VariableNames',{'IdE','Mw','Kinematics'});
                eventInfo = [eventInfo;T]; %#ok<AGROW>
            end
        else
            warning('Event-list file not found: %s',filename);
        end
    end

    if ~isempty(eventInfo)
        [~,ia] = unique(eventInfo.IdE,'stable');
        eventInfo = eventInfo(ia,:);
    end
end

function mask = matchAttributeValue(S,fieldName,targetValue)

    raw = {S.(fieldName)}';
    mask = false(numel(raw),1);
    targetText = strtrim(string(targetValue));
    targetNumeric = str2double(targetText);
    targetIsNumeric = isfinite(targetNumeric);

    for i = 1:numel(raw)
        value = raw{i};

        if isempty(value)
            continue;
        end

        valueText = strtrim(string(value));

        % Text comparison is case-insensitive. This handles values such as
        % 'Main', 'MAIN', 'Principal', etc.
        if strcmpi(valueText,targetText)
            mask(i) = true;
            continue;
        end

        % If both values are numerically interpretable, compare numerically
        % so that 1, 1.0 and '1' are treated as equivalent.
        if targetIsNumeric
            valueNumeric = str2double(valueText);
            if isfinite(valueNumeric) && abs(valueNumeric-targetNumeric)<=1e-9
                mask(i) = true;
            end
        end
    end
end

function values = getNumericField(S,fieldName)

    raw = {S.(fieldName)}';
    values = nan(numel(raw),1);

    for i = 1:numel(raw)
        if isnumeric(raw{i})
            values(i) = double(raw{i});
        else
            values(i) = str2double(string(raw{i}));
        end
    end
end

function [features,stats] = cleanObservedRank1Geometry( ...
    features,duplicateTolerance,minimumLength)

    % Geometry-preserving cleaning only. No endpoint snapping, feature
    % merging, coordinate averaging, or Douglas-Peucker replacement occurs.
    stats = struct( ...
        'RemovedVertices',0, ...
        'SnappedTips',0, ...
        'MergedPairs',0);

    if isempty(features)
        return;
    end

    cleaned = repmat(features(1),0,1);

    for i = 1:numel(features)

        [x,y,nRemoved] = cleanPolylineVertices( ...
            features(i).x,features(i).y,duplicateTolerance);

        stats.RemovedVertices = stats.RemovedVertices+nRemoved;

        if numel(x)<2
            continue;
        end

        if sum(hypot(diff(x),diff(y)))<minimumLength
            continue;
        end

        f = features(i);
        f.x = x;
        f.y = y;
        f = refreshFeatureGeometry(f);
        cleaned(end+1,1) = f; %#ok<AGROW>
    end

    features = cleaned;
end

function [features,stats] = preprocessRank1Geometry( ...
    features,duplicateTolerance,tipSnapTolerance,mergeTolerance, ...
    maximumMergeTurn,maximumTortuosity,rejectSelfIntersections, ...
    maximumMergePasses,minimumLength)

    stats = struct( ...
        'RemovedVertices',0, ...
        'SnappedTips',0, ...
        'MergedPairs',0);

    % 1. Remove duplicate and zero-length consecutive vertices.
    cleaned = repmat(features(1),0,1);

    for i = 1:numel(features)

        [x,y,nRemoved] = cleanPolylineVertices( ...
            features(i).x,features(i).y,duplicateTolerance);

        stats.RemovedVertices = stats.RemovedVertices+nRemoved;

        if numel(x)<2
            continue;
        end

        L = sum(hypot(diff(x),diff(y)));

        if L<minimumLength
            continue;
        end

        f = features(i);
        f.x = x;
        f.y = y;
        f = refreshFeatureGeometry(f);
        cleaned(end+1,1) = f; %#ok<AGROW>
    end

    features = cleaned;

    if isempty(features)
        return;
    end

    % 2. Snap nearby original tips to common cluster centres.
    [features,nSnapped] = snapFeatureTips(features,tipSnapTolerance);
    stats.SnappedTips = nSnapped;

    % 3. Merge nearly continuous features. Repeat because one merge can
    % create a new admissible endpoint pair.
    for iPass = 1:maximumMergePasses

        [found,iFeature,jFeature,endI,endJ] = ...
            findBestMergePair( ...
                features,mergeTolerance,maximumMergeTurn, ...
                maximumTortuosity,rejectSelfIntersections);

        if ~found
            break;
        end

        mergedFeature = mergeTwoFeatures( ...
            features(iFeature),features(jFeature),endI,endJ);

        keep = true(numel(features),1);
        keep([iFeature jFeature]) = false;
        features = features(keep);
        features(end+1,1) = mergedFeature; %#ok<AGROW>
        stats.MergedPairs = stats.MergedPairs+1;

        % Re-snap after each merge to avoid accumulated sub-metre offsets.
        [features,nAdditionalSnap] = ...
            snapFeatureTips(features,tipSnapTolerance);
        stats.SnappedTips = stats.SnappedTips+nAdditionalSnap;
    end

    % 4. Final cleanup and geometry refresh.
    finalFeatures = repmat(features(1),0,1);

    for i = 1:numel(features)

        [x,y,nRemoved] = cleanPolylineVertices( ...
            features(i).x,features(i).y,duplicateTolerance);

        stats.RemovedVertices = stats.RemovedVertices+nRemoved;

        if numel(x)<2 || sum(hypot(diff(x),diff(y)))<minimumLength
            continue;
        end

        f = features(i);
        f.x = x;
        f.y = y;
        f = refreshFeatureGeometry(f);
        finalFeatures(end+1,1) = f; %#ok<AGROW>
    end

    features = finalFeatures;
end

function [x,y,nRemoved] = cleanPolylineVertices(x,y,tolerance)

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x)<2
        nRemoved = 0;
        return;
    end

    keep = true(numel(x),1);
    lastKept = 1;

    for i = 2:numel(x)
        if hypot(x(i)-x(lastKept),y(i)-y(lastKept))<=tolerance
            keep(i) = false;
        else
            lastKept = i;
        end
    end

    nRemoved = sum(~keep);
    x = x(keep);
    y = y(keep);
end

function feature = refreshFeatureGeometry(feature)

    feature.x = feature.x(:);
    feature.y = feature.y(:);
    feature.cum_m = polylineCumulativeDistance(feature.x,feature.y);
    feature.length_m = feature.cum_m(end);
    feature.startXY = [feature.x(1),feature.y(1)];
    feature.endXY = [feature.x(end),feature.y(end)];
end

function [features,nSnapped] = snapFeatureTips(features,tolerance)

    nFeatures = numel(features);
    points = zeros(2*nFeatures,2);

    for i = 1:nFeatures
        points(2*i-1,:) = [features(i).x(1),features(i).y(1)];
        points(2*i,:)   = [features(i).x(end),features(i).y(end)];
    end

    [centres,labels] = clusterPoints(points,tolerance);
    nSnapped = 0;

    for i = 1:nFeatures

        startCentre = centres(labels(2*i-1),:);
        endCentre   = centres(labels(2*i),:);

        if hypot(features(i).x(1)-startCentre(1), ...
                 features(i).y(1)-startCentre(2))>1e-9
            nSnapped = nSnapped+1;
        end

        if hypot(features(i).x(end)-endCentre(1), ...
                 features(i).y(end)-endCentre(2))>1e-9
            nSnapped = nSnapped+1;
        end

        features(i).x(1) = startCentre(1);
        features(i).y(1) = startCentre(2);
        features(i).x(end) = endCentre(1);
        features(i).y(end) = endCentre(2);
        features(i) = refreshFeatureGeometry(features(i));
    end
end

function [found,bestI,bestJ,bestEndI,bestEndJ] = ...
    findBestMergePair(features,distanceTolerance,maximumTurn, ...
    maximumTortuosity,rejectSelfIntersections)

    found = false;
    bestI = NaN;
    bestJ = NaN;
    bestEndI = NaN;
    bestEndJ = NaN;
    bestDistance = Inf;
    bestTurn = Inf;

    for i = 1:numel(features)-1
        for j = i+1:numel(features)
            for endI = 1:2
                for endJ = 1:2

                    pointI = featureEndpoint(features(i),endI);
                    pointJ = featureEndpoint(features(j),endJ);
                    distance = hypot(pointI(1)-pointJ(1),pointI(2)-pointJ(2));

                    if distance>distanceTolerance
                        continue;
                    end

                    approachI = directionTowardEndpoint(features(i),endI);
                    departJ = directionAwayFromEndpoint(features(j),endJ);
                    turn = vectorAngleDegrees(approachI,departJ);

                    if turn>maximumTurn
                        continue;
                    end

                    proposedMerge = mergeTwoFeatures( ...
                        features(i),features(j),endI,endJ);

                    endpointDistance = hypot( ...
                        proposedMerge.x(end)-proposedMerge.x(1), ...
                        proposedMerge.y(end)-proposedMerge.y(1));

                    mergeTortuosity = ...
                        proposedMerge.length_m/max(endpointDistance,eps);

                    if mergeTortuosity>maximumTortuosity
                        continue;
                    end

                    if rejectSelfIntersections && ...
                            polylineHasSelfIntersection( ...
                                proposedMerge.x,proposedMerge.y)
                        continue;
                    end

                    if distance<bestDistance-1e-9 || ...
                            (abs(distance-bestDistance)<=1e-9 && turn<bestTurn)

                        found = true;
                        bestI = i;
                        bestJ = j;
                        bestEndI = endI;
                        bestEndJ = endJ;
                        bestDistance = distance;
                        bestTurn = turn;
                    end
                end
            end
        end
    end
end

function point = featureEndpoint(feature,endCode)

    if endCode==1
        point = [feature.x(1),feature.y(1)];
    else
        point = [feature.x(end),feature.y(end)];
    end
end

function direction = directionTowardEndpoint(feature,endCode)

    if endCode==1
        direction = [ ...
            feature.x(1)-feature.x(2), ...
            feature.y(1)-feature.y(2)];
    else
        direction = [ ...
            feature.x(end)-feature.x(end-1), ...
            feature.y(end)-feature.y(end-1)];
    end

    direction = normalizeVectorLocal(direction);
end

function direction = directionAwayFromEndpoint(feature,endCode)

    if endCode==1
        direction = [ ...
            feature.x(2)-feature.x(1), ...
            feature.y(2)-feature.y(1)];
    else
        direction = [ ...
            feature.x(end-1)-feature.x(end), ...
            feature.y(end-1)-feature.y(end)];
    end

    direction = normalizeVectorLocal(direction);
end

function merged = mergeTwoFeatures(featureA,featureB,endA,endB)

    % Orient A so the selected endpoint is last.
    if endA==1
        xA = flipud(featureA.x(:));
        yA = flipud(featureA.y(:));
    else
        xA = featureA.x(:);
        yA = featureA.y(:);
    end

    % Orient B so the selected endpoint is first.
    if endB==2
        xB = flipud(featureB.x(:));
        yB = flipud(featureB.y(:));
    else
        xB = featureB.x(:);
        yB = featureB.y(:);
    end

    junction = 0.5*([xA(end),yA(end)]+[xB(1),yB(1)]);
    xA(end) = junction(1);
    yA(end) = junction(2);
    xB(1) = junction(1);
    yB(1) = junction(2);

    merged = featureA;
    merged.x = [xA;xB(2:end)];
    merged.y = [yA;yB(2:end)];

    if isequaln(featureA.IdS,featureB.IdS)
        merged.IdS = featureA.IdS;
    else
        merged.IdS = -1;
    end

    merged.FeatureIndex = min(featureA.FeatureIndex,featureB.FeatureIndex);
    merged.PartWithinFeature = 1;
    merged = refreshFeatureGeometry(merged);
end

function tf = polylineHasSelfIntersection(x,y)

    x = x(:);
    y = y(:);
    nSegments = numel(x)-1;
    tf = false;

    for i = 1:nSegments-2
        p1 = [x(i) y(i)];
        p2 = [x(i+1) y(i+1)];

        for j = i+2:nSegments

            % Adjacent segments share a legal vertex. Also skip the first
            % and last segment when they only meet at a closed endpoint.
            if j==i+1
                continue;
            end

            q1 = [x(j) y(j)];
            q2 = [x(j+1) y(j+1)];

            if segmentsIntersectProperly(p1,p2,q1,q2)
                tf = true;
                return;
            end
        end
    end
end

function tf = segmentsIntersectProperly(p1,p2,q1,q2)

    tolerance = 1e-9;

    o1 = cross2d(p2-p1,q1-p1);
    o2 = cross2d(p2-p1,q2-p1);
    o3 = cross2d(q2-q1,p1-q1);
    o4 = cross2d(q2-q1,p2-q1);

    tf = ...
        ((o1>tolerance && o2<-tolerance) || ...
         (o1<-tolerance && o2>tolerance)) && ...
        ((o3>tolerance && o4<-tolerance) || ...
         (o3<-tolerance && o4>tolerance));
end

function value = cross2d(a,b)
    value = a(1)*b(2)-a(2)*b(1);
end

function angle = vectorAngleDegrees(v1,v2)

    v1 = normalizeVectorLocal(v1);
    v2 = normalizeVectorLocal(v2);
    cosine = max(-1,min(1,dot(v1,v2)));
    angle = acosd(cosine);
end

function vector = normalizeVectorLocal(vector)

    normValue = hypot(vector(1),vector(2));

    if normValue<=eps
        vector = [1 0];
    else
        vector = vector/normValue;
    end
end

function parts = extractPolylineParts(PF)

    parts = struct('lon',{},'lat',{},'IdS',{}, ...
        'FeatureIndex',{},'PartWithinFeature',{});

    counter = 0;

    for i = 1:numel(PF)

        lon = PF(i).X(:);
        lat = PF(i).Y(:);
        valid = ~(isnan(lon) | isnan(lat));
        validIdx = find(valid);

        if isempty(validIdx)
            continue;
        end

        breaks = find(diff(validIdx)>1);
        starts = [validIdx(1);validIdx(breaks+1)];
        ends   = [validIdx(breaks);validIdx(end)];

        for j = 1:numel(starts)

            idx = starts(j):ends(j);

            if numel(idx)<2
                continue;
            end

            counter = counter+1;
            parts(counter).lon = lon(idx);
            parts(counter).lat = lat(idx);

            if isfield(PF,'IdS')
                if isnumeric(PF(i).IdS)
                    parts(counter).IdS = double(PF(i).IdS);
                else
                    parts(counter).IdS = str2double(string(PF(i).IdS));
                end
            else
                parts(counter).IdS = i;
            end

            parts(counter).FeatureIndex = i;
            parts(counter).PartWithinFeature = j;
        end
    end
end

function [x,y] = lonlat2local(lon,lat,lon0,lat0)

    R = 6371008.8;
    x = deg2rad(lon-lon0).*R.*cosd(lat0);
    y = deg2rad(lat-lat0).*R;
end

function [lon,lat] = local2lonlat(x,y,lon0,lat0)

    R = 6371008.8;
    lon = lon0 + rad2deg(x./(R*cosd(lat0)));
    lat = lat0 + rad2deg(y./R);
end

function cum = polylineCumulativeDistance(x,y)

    cum = [0;cumsum(hypot(diff(x),diff(y)))];
end

function [splitPositions,junctionTable] = findDynamicJunctions( ...
    original,junctionTolerance,minimumSplitDistance,mergeTolerance)

    n = numel(original);
    splitPositions = cell(n,1);

    for i = 1:n
        splitPositions{i} = [0;original(i).length_m];
    end

    TipPart = [];
    TipEnd = [];
    TipX = [];
    TipY = [];
    ReceiverPart = [];
    AlongReceiver_m = [];
    Distance_m = [];
    ProjectionX = [];
    ProjectionY = [];

    for iTipPart = 1:n

        tipCoordinates = [ ...
            original(iTipPart).startXY; ...
            original(iTipPart).endXY];

        for iEnd = 1:2

            qx = tipCoordinates(iEnd,1);
            qy = tipCoordinates(iEnd,2);

            for j = 1:n

                if j == iTipPart
                    continue;
                end

                [d,s,xp,yp] = projectPointOnPolyline( ...
                    qx,qy,original(j).x,original(j).y,original(j).cum_m);

                if d > junctionTolerance
                    continue;
                end

                if s <= minimumSplitDistance || ...
                        s >= original(j).length_m-minimumSplitDistance
                    continue;
                end

                splitPositions{j}(end+1,1) = s; %#ok<AGROW>

                TipPart(end+1,1) = iTipPart; %#ok<AGROW>
                TipEnd(end+1,1) = iEnd; %#ok<AGROW>
                TipX(end+1,1) = qx; %#ok<AGROW>
                TipY(end+1,1) = qy; %#ok<AGROW>
                ReceiverPart(end+1,1) = j; %#ok<AGROW>
                AlongReceiver_m(end+1,1) = s; %#ok<AGROW>
                Distance_m(end+1,1) = d; %#ok<AGROW>
                ProjectionX(end+1,1) = xp; %#ok<AGROW>
                ProjectionY(end+1,1) = yp; %#ok<AGROW>
            end
        end
    end

    for i = 1:n
        splitPositions{i} = mergeCloseValues( ...
            sort(splitPositions{i}),mergeTolerance);
    end

    junctionTable = table(TipPart,TipEnd,TipX,TipY,ReceiverPart, ...
        AlongReceiver_m,Distance_m,ProjectionX,ProjectionY);
end

function values = mergeCloseValues(values,tolerance)

    if isempty(values)
        return;
    end

    values = sort(values(:));
    merged = values(1);

    for i = 2:numel(values)
        if values(i)-merged(end) <= tolerance
            merged(end) = mean([merged(end),values(i)]);
        else
            merged(end+1,1) = values(i); %#ok<AGROW>
        end
    end

    values = merged;
end

function [splitPositions,records] = addDouglasPeuckerSplitPositions( ...
    original,splitPositions,mergeTolerance,minTolerance,maxTolerance,lengthFraction)

    OriginalPartIndex = zeros(0,1);
    VertexIndex = zeros(0,1);
    AlongFeature_m = zeros(0,1);
    X_m = zeros(0,1);
    Y_m = zeros(0,1);
    Tolerance_m = zeros(0,1);

    for i = 1:numel(original)

        tolerance = max(minTolerance, ...
            min(maxTolerance,lengthFraction*original(i).length_m));

        indices = douglasPeuckerIndices( ...
            original(i).x,original(i).y,tolerance);

        % Tips already exist. Only internal retained vertices are added.
        indices = indices(indices>1 & indices<numel(original(i).x));

        for k = 1:numel(indices)
            iv = indices(k);
            along = original(i).cum_m(iv);

            splitPositions{i}(end+1,1) = along; %#ok<AGROW>

            OriginalPartIndex(end+1,1) = i; %#ok<AGROW>
            VertexIndex(end+1,1) = iv; %#ok<AGROW>
            AlongFeature_m(end+1,1) = along; %#ok<AGROW>
            X_m(end+1,1) = original(i).x(iv); %#ok<AGROW>
            Y_m(end+1,1) = original(i).y(iv); %#ok<AGROW>
            Tolerance_m(end+1,1) = tolerance; %#ok<AGROW>
        end

        splitPositions{i} = mergeCloseValues( ...
            sort(splitPositions{i}),mergeTolerance);
    end

    records = table(OriginalPartIndex,VertexIndex,AlongFeature_m, ...
        X_m,Y_m,Tolerance_m);
end

function mask = identifySignificantExitNodes(nodes,records,tolerance)

    mask = false(height(nodes),1);

    for i = 1:height(records)
        d = hypot(nodes.X_m-records.X_m(i), ...
                  nodes.Y_m-records.Y_m(i));
        [dmin,j] = min(d);
        if dmin<=tolerance
            mask(j) = true;
        end
    end
end

function indices = douglasPeuckerIndices(x,y,tolerance)

    x = x(:);
    y = y(:);
    n = numel(x);

    if n<=2 || tolerance<=0
        indices = (1:n)';
        return;
    end

    keep = false(n,1);
    keep([1 n]) = true;

    stackStart = 1;
    stackEnd = n;

    while ~isempty(stackStart)

        i1 = stackStart(end);
        i2 = stackEnd(end);
        stackStart(end) = [];
        stackEnd(end) = [];

        if i2<=i1+1
            continue;
        end

        distance = pointToSegmentDistanceDP( ...
            x(i1+1:i2-1),y(i1+1:i2-1), ...
            x(i1),y(i1),x(i2),y(i2));

        [maximumDistance,localIndex] = max(distance);

        if maximumDistance>tolerance
            splitIndex = i1+localIndex;
            keep(splitIndex) = true;

            stackStart(end+1,1) = i1; %#ok<AGROW>
            stackEnd(end+1,1) = splitIndex; %#ok<AGROW>
            stackStart(end+1,1) = splitIndex; %#ok<AGROW>
            stackEnd(end+1,1) = i2; %#ok<AGROW>
        end
    end

    indices = find(keep);
end

function distance = pointToSegmentDistanceDP(px,py,x1,y1,x2,y2)

    vx = x2-x1;
    vy = y2-y1;
    lengthSquared = vx^2+vy^2;

    if lengthSquared<=0
        distance = hypot(px-x1,py-y1);
        return;
    end

    t = ((px-x1)*vx+(py-y1)*vy)/lengthSquared;
    t = max(0,min(1,t));

    projectionX = x1+t*vx;
    projectionY = y1+t*vy;

    distance = hypot(px-projectionX,py-projectionY);
end

function edges = splitOriginalFeatures(original,splitPositions,tolerance)

    edges = struct('OriginalPartIndex',{},'IdS',{},'x',{},'y',{}, ...
        'StartAlong_m',{},'EndAlong_m',{},'Length_m',{});

    counter = 0;

    for i = 1:numel(original)

        cuts = splitPositions{i};
        cuts(1) = 0;
        cuts(end) = original(i).length_m;
        cuts = mergeCloseValues(cuts,tolerance);

        for k = 1:numel(cuts)-1

            s1 = cuts(k);
            s2 = cuts(k+1);

            if s2-s1 <= tolerance
                continue;
            end

            [xs,ys] = slicePolylineByDistance( ...
                original(i).x,original(i).y,original(i).cum_m,s1,s2);

            if numel(xs)<2
                continue;
            end

            counter = counter+1;
            edges(counter).OriginalPartIndex = i;
            edges(counter).IdS = original(i).IdS;
            edges(counter).x = xs;
            edges(counter).y = ys;
            edges(counter).StartAlong_m = s1;
            edges(counter).EndAlong_m = s2;
            edges(counter).Length_m = s2-s1;
        end
    end
end

function [xs,ys] = slicePolylineByDistance(x,y,cum,s1,s2)

    [x1,y1] = pointAtDistance(x,y,cum,s1);
    [x2,y2] = pointAtDistance(x,y,cum,s2);

    internal = find(cum>s1 & cum<s2);

    xs = [x1;x(internal);x2];
    ys = [y1;y(internal);y2];

    keep = [true;hypot(diff(xs),diff(ys))>1e-9];
    xs = xs(keep);
    ys = ys(keep);
end

function [xp,yp] = pointAtDistance(x,y,cum,s)

    if s <= 0
        xp = x(1); yp = y(1); return;
    end

    if s >= cum(end)
        xp = x(end); yp = y(end); return;
    end

    k = find(cum<=s,1,'last');

    if k == numel(cum)
        xp = x(end); yp = y(end); return;
    end

    ds = cum(k+1)-cum(k);

    if ds <= 0
        t = 0;
    else
        t = (s-cum(k))/ds;
    end

    xp = x(k)+t*(x(k+1)-x(k));
    yp = y(k)+t*(y(k+1)-y(k));
end

function [distance,along,xp,yp] = projectPointOnPolyline(qx,qy,x,y,cum)

    bestD2 = Inf;
    along = NaN;
    xp = NaN;
    yp = NaN;

    for k = 1:numel(x)-1

        vx = x(k+1)-x(k);
        vy = y(k+1)-y(k);
        len2 = vx^2+vy^2;

        if len2 <= 0
            continue;
        end

        t = ((qx-x(k))*vx+(qy-y(k))*vy)/len2;
        t = max(0,min(1,t));

        xn = x(k)+t*vx;
        yn = y(k)+t*vy;
        d2 = (qx-xn)^2+(qy-yn)^2;

        if d2<bestD2
            bestD2 = d2;
            xp = xn;
            yp = yn;
            along = cum(k)+t*sqrt(len2);
        end
    end

    distance = sqrt(bestD2);
end

function [nodes,traceEdges] = buildSnappedTraceGraph(rawEdges,snapTol,original)

    endpointXY = zeros(2*numel(rawEdges),2);

    for i = 1:numel(rawEdges)
        endpointXY(2*i-1,:) = [rawEdges(i).x(1),rawEdges(i).y(1)];
        endpointXY(2*i,:)   = [rawEdges(i).x(end),rawEdges(i).y(end)];
    end

    [nodeXY,endpointNode] = clusterPoints(endpointXY,snapTol);

    nNodes = size(nodeXY,1);
    nodes = table((1:nNodes)',nodeXY(:,1),nodeXY(:,2), ...
        'VariableNames',{'NodeID','X_m','Y_m'});

    nE = numel(rawEdges);

    EdgeID = (1:nE)';
    Node1 = zeros(nE,1);
    Node2 = zeros(nE,1);
    EdgeType = repmat("TRACE",nE,1);
    Length_m = zeros(nE,1);
    OriginalPartIndex = zeros(nE,1);
    IdS = zeros(nE,1);
    StartAlong_m = zeros(nE,1);
    EndAlong_m = zeros(nE,1);
    XGeom = cell(nE,1);
    YGeom = cell(nE,1);

    for i = 1:nE
        Node1(i) = endpointNode(2*i-1);
        Node2(i) = endpointNode(2*i);
        Length_m(i) = rawEdges(i).Length_m;
        OriginalPartIndex(i) = rawEdges(i).OriginalPartIndex;
        IdS(i) = rawEdges(i).IdS;
        StartAlong_m(i) = rawEdges(i).StartAlong_m;
        EndAlong_m(i) = rawEdges(i).EndAlong_m;
        XGeom{i} = rawEdges(i).x(:);
        YGeom{i} = rawEdges(i).y(:);
    end

    traceEdges = table(EdgeID,Node1,Node2,EdgeType,Length_m, ...
        OriginalPartIndex,IdS,StartAlong_m,EndAlong_m,XGeom,YGeom);

    % Mark nodes that coincide with original feature tips.
    nodes.IsOriginalTip = false(height(nodes),1);

    originalTips = [];
    for i = 1:numel(original)
        originalTips = [originalTips;original(i).startXY;original(i).endXY]; %#ok<AGROW>
    end

    for i = 1:size(originalTips,1)
        d = hypot(nodes.X_m-originalTips(i,1),nodes.Y_m-originalTips(i,2));
        [dmin,j] = min(d);
        if dmin<=snapTol
            nodes.IsOriginalTip(j) = true;
        end
    end
end

function [centres,labels] = clusterPoints(points,tolerance)

    n = size(points,1);
    labels = zeros(n,1);
    centres = zeros(0,2);
    nClusters = 0;

    for i = 1:n

        if nClusters == 0
            nClusters = 1;
            centres(1,:) = points(i,:);
            labels(i) = 1;
            continue;
        end

        d = hypot(centres(:,1)-points(i,1),centres(:,2)-points(i,2));
        [dmin,j] = min(d);

        if dmin<=tolerance
            labels(i) = j;
            members = points(labels==j,:);
            centres(j,:) = mean(members,1);
        else
            nClusters = nClusters+1;
            centres(nClusters,:) = points(i,:); %#ok<AGROW>
            labels(i) = nClusters;
        end
    end
end

function junctionEdges = buildJunctionConnectorEdges( ...
    nodes,junctionRecords,snapTolerance)

    Node1 = zeros(0,1);
    Node2 = zeros(0,1);
    Length_m = zeros(0,1);
    XGeom = cell(0,1);
    YGeom = cell(0,1);

    for i = 1:height(junctionRecords)

        dTip = hypot(nodes.X_m-junctionRecords.TipX(i), ...
                     nodes.Y_m-junctionRecords.TipY(i));
        [tipDistance,tipNode] = min(dTip);

        dProjection = hypot(nodes.X_m-junctionRecords.ProjectionX(i), ...
                            nodes.Y_m-junctionRecords.ProjectionY(i));
        [projectionDistance,projectionNode] = min(dProjection);

        % The split and original-tip nodes should normally be very close to
        % the recorded coordinates. A slightly relaxed check is used to
        % accommodate node averaging during snapping.
        if tipDistance>2*snapTolerance || projectionDistance>2*snapTolerance
            continue;
        end

        if tipNode==projectionNode
            continue;
        end

        pair = sort([tipNode projectionNode]);

        if ~isempty(Node1)
            existing = sort([Node1 Node2],2);
            if any(existing(:,1)==pair(1) & existing(:,2)==pair(2))
                continue;
            end
        end

        Node1(end+1,1) = tipNode; %#ok<AGROW>
        Node2(end+1,1) = projectionNode; %#ok<AGROW>
        Length_m(end+1,1) = hypot( ...
            nodes.X_m(tipNode)-nodes.X_m(projectionNode), ...
            nodes.Y_m(tipNode)-nodes.Y_m(projectionNode)); %#ok<AGROW>
        XGeom{end+1,1} = [nodes.X_m(tipNode);nodes.X_m(projectionNode)]; %#ok<AGROW>
        YGeom{end+1,1} = [nodes.Y_m(tipNode);nodes.Y_m(projectionNode)]; %#ok<AGROW>
    end

    n = numel(Node1);
    EdgeID = (1:n)';
    EdgeType = repmat("JUNCTION",n,1);
    OriginalPartIndex = nan(n,1);
    IdS = nan(n,1);
    StartAlong_m = nan(n,1);
    EndAlong_m = nan(n,1);

    junctionEdges = table(EdgeID,Node1,Node2,EdgeType,Length_m, ...
        OriginalPartIndex,IdS,StartAlong_m,EndAlong_m,XGeom,YGeom);
end

function gapEdges = buildArtificialGapEdges( ...
    nodes,existingEdges,maxGap,minGap,mode,maxConnectionsPerExit)

    existingPairs = sort([existingEdges.Node1 existingEdges.Node2],2);

    Node1 = zeros(0,1);
    Node2 = zeros(0,1);
    Length_m = zeros(0,1);
    XGeom = cell(0,1);
    YGeom = cell(0,1);

    if nargin<6 || isempty(maxConnectionsPerExit)
        maxConnectionsPerExit = Inf;
    end

    % Original-part membership of every node. This prevents artificial
    % shortcuts from an internal DP vertex to a tip of the same feature.
    nodeParts = cell(height(nodes),1);

    traceMask = existingEdges.EdgeType=="TRACE";
    traceRows = find(traceMask);

    for ii = 1:numel(traceRows)
        e = traceRows(ii);
        p = existingEdges.OriginalPartIndex(e);
        if isfinite(p)
            nodeParts{existingEdges.Node1(e)}(end+1) = p; %#ok<AGROW>
            nodeParts{existingEdges.Node2(e)}(end+1) = p; %#ok<AGROW>
        end
    end

    for i = 1:height(nodes)
        nodeParts{i} = unique(nodeParts{i});
    end

    switch lower(mode)
        case 'significanttoterminal'
            startNodes = find(nodes.IsSignificantExit);

            for is = 1:numel(startNodes)
                i = startNodes(is);

                candidateTips = find(nodes.IsOriginalTip);
                candidateTips(candidateTips==i) = [];

                if isempty(candidateTips)
                    continue;
                end

                d = hypot(nodes.X_m(candidateTips)-nodes.X_m(i), ...
                          nodes.Y_m(candidateTips)-nodes.Y_m(i));

                valid = d>=minGap & d<=maxGap;

                % Avoid connection within the same original feature.
                for k = 1:numel(candidateTips)
                    if ~isempty(intersect(nodeParts{i},nodeParts{candidateTips(k)}))
                        valid(k) = false;
                    end
                end

                candidateTips = candidateTips(valid);
                d = d(valid);

                if isempty(candidateTips)
                    continue;
                end

                [d,order] = sort(d,'ascend');
                candidateTips = candidateTips(order);

                nKeep = min(maxConnectionsPerExit,numel(candidateTips));
                candidateTips = candidateTips(1:nKeep);
                d = d(1:nKeep);

                for k = 1:numel(candidateTips)
                    j = candidateTips(k);
                    pair = sort([i j]);

                    if any(existingPairs(:,1)==pair(1) & ...
                           existingPairs(:,2)==pair(2))
                        continue;
                    end

                    if ~isempty(Node1)
                        createdPairs = sort([Node1 Node2],2);
                        if any(createdPairs(:,1)==pair(1) & ...
                               createdPairs(:,2)==pair(2))
                            continue;
                        end
                    end

                    Node1(end+1,1) = i; %#ok<AGROW>
                    Node2(end+1,1) = j; %#ok<AGROW>
                    Length_m(end+1,1) = d(k); %#ok<AGROW>
                    XGeom{end+1,1} = [nodes.X_m(i);nodes.X_m(j)]; %#ok<AGROW>
                    YGeom{end+1,1} = [nodes.Y_m(i);nodes.Y_m(j)]; %#ok<AGROW>
                end
            end

        otherwise
            for i = 1:height(nodes)-1
                for j = i+1:height(nodes)

                    switch lower(mode)
                        case 'terminalonly'
                            allowed = nodes.IsTerminal(i) && nodes.IsTerminal(j);
                        case 'terminaltoall'
                            allowed = nodes.IsTerminal(i) || nodes.IsTerminal(j);
                        case 'all'
                            allowed = true;
                        otherwise
                            error('Unknown gapConnectionMode: %s',mode);
                    end

                    if ~allowed
                        continue;
                    end

                    pair = [i j];
                    if any(existingPairs(:,1)==pair(1) & ...
                           existingPairs(:,2)==pair(2))
                        continue;
                    end

                    d = hypot(nodes.X_m(i)-nodes.X_m(j), ...
                              nodes.Y_m(i)-nodes.Y_m(j));

                    if d<minGap || d>maxGap
                        continue;
                    end

                    Node1(end+1,1) = i; %#ok<AGROW>
                    Node2(end+1,1) = j; %#ok<AGROW>
                    Length_m(end+1,1) = d; %#ok<AGROW>
                    XGeom{end+1,1} = [nodes.X_m(i);nodes.X_m(j)]; %#ok<AGROW>
                    YGeom{end+1,1} = [nodes.Y_m(i);nodes.Y_m(j)]; %#ok<AGROW>
                end
            end
    end

    n = numel(Node1);
    EdgeID = (1:n)';
    EdgeType = repmat("GAP",n,1);
    OriginalPartIndex = nan(n,1);
    IdS = nan(n,1);
    StartAlong_m = nan(n,1);
    EndAlong_m = nan(n,1);

    gapEdges = table(EdgeID,Node1,Node2,EdgeType,Length_m, ...
        OriginalPartIndex,IdS,StartAlong_m,EndAlong_m,XGeom,YGeom);
end

function allEdges = concatenateEdgeTables(traceEdges,gapEdges)

    allEdges = [traceEdges;gapEdges];
    allEdges.EdgeID = (1:height(allEdges))';
end

function tipNodes = identifyOriginalTipNodes(original,nodes,tolerance)

    TipID = [];
    NodeID = [];
    OriginalPartIndex = [];
    EndCode = [];
    X_m = [];
    Y_m = [];

    counter = 0;

    for i = 1:numel(original)

        coordinates = [original(i).startXY;original(i).endXY];

        for e = 1:2

            d = hypot(nodes.X_m-coordinates(e,1),nodes.Y_m-coordinates(e,2));
            [dmin,j] = min(d);

            if dmin>tolerance
                continue;
            end

            counter = counter+1;
            TipID(counter,1) = counter; %#ok<AGROW>
            NodeID(counter,1) = j; %#ok<AGROW>
            OriginalPartIndex(counter,1) = i; %#ok<AGROW>
            EndCode(counter,1) = e; %#ok<AGROW>
            X_m(counter,1) = nodes.X_m(j); %#ok<AGROW>
            Y_m(counter,1) = nodes.Y_m(j); %#ok<AGROW>
        end
    end

    tipNodes = table(TipID,NodeID,OriginalPartIndex,EndCode,X_m,Y_m);
    [~,ia] = unique(tipNodes.NodeID,'stable');
    tipNodes = tipNodes(ia,:);
end

function [candidateTips,pairs,startTips,endTips] = ...
    chooseOppositeExtremeTipPairs( ...
        tipNodes,globalNodeProgress,nExtremePerSide,maxPairs)

    % Keep one row per graph node.
    [~,ia] = unique(tipNodes.NodeID,'stable');
    uniqueTips = tipNodes(ia,:);

    if isempty(uniqueTips)
        candidateTips = uniqueTips;
        startTips = uniqueTips;
        endTips = uniqueTips;
        pairs = table( ...
            zeros(0,1),zeros(0,1),zeros(0,1), ...
            'VariableNames',{'StartNode','EndNode','ProgressSpan_m'});
        return;
    end

    progress = globalNodeProgress(uniqueTips.NodeID);
    [~,order] = sort(progress,'ascend');
    uniqueTips = uniqueTips(order,:);
    progress = progress(order);

    nAvailable = height(uniqueTips);
    nSide = min(nExtremePerSide,floor(nAvailable/2));

    if nSide<1
        candidateTips = uniqueTips;
        startTips = uniqueTips([],:);
        endTips = uniqueTips([],:);
        pairs = table( ...
            zeros(0,1),zeros(0,1),zeros(0,1), ...
            'VariableNames',{'StartNode','EndNode','ProgressSpan_m'});
        return;
    end

    startTips = uniqueTips(1:nSide,:);
    endTips = uniqueTips(end-nSide+1:end,:);

    % In unusual compact geometries, the two pools can share nodes.
    % Remove shared nodes from the end-side pool so that every pair joins
    % genuinely different graph nodes.
    sharedNodeIDs = intersect(startTips.NodeID,endTips.NodeID);
    if ~isempty(sharedNodeIDs)
        endTips = endTips(~ismember(endTips.NodeID,sharedNodeIDs),:);
    end

    candidateTips = [startTips;endTips];
    [~,ia] = unique(candidateTips.NodeID,'stable');
    candidateTips = candidateTips(ia,:);

    if isempty(startTips) || isempty(endTips)
        pairs = table( ...
            zeros(0,1),zeros(0,1),zeros(0,1), ...
            'VariableNames',{'StartNode','EndNode','ProgressSpan_m'});
        return;
    end

    [startGrid,endGrid] = ndgrid(startTips.NodeID,endTips.NodeID);

    StartNode = startGrid(:);
    EndNode = endGrid(:);

    valid = StartNode~=EndNode;
    StartNode = StartNode(valid);
    EndNode = EndNode(valid);

    ProgressSpan_m = abs( ...
        globalNodeProgress(EndNode)- ...
        globalNodeProgress(StartNode));

    pairs = table(StartNode,EndNode,ProgressSpan_m);
    pairs = unique(pairs,'rows','stable');
    pairs = sortrows(pairs,'ProgressSpan_m','descend');

    if isfinite(maxPairs) && height(pairs)>maxPairs
        pairs = pairs(1:maxPairs,:);
    end
end

function candidate = solveDynamicStatePath( ...
    nodes,edges,startNode,endNode,W,preferredGap,maxTurn)

    nEdges = height(edges);
    nStates = 2*nEdges;

    stateEdge = repelem((1:nEdges)',2);
    stateDirection = repmat([1;-1],nEdges,1);
    stateStartNode = zeros(nStates,1);
    stateEndNode = zeros(nStates,1);
    stateStartAz = zeros(nStates,1);
    stateEndAz = zeros(nStates,1);
    stateTraversalCost = zeros(nStates,1);

    globalVec = [nodes.X_m(endNode)-nodes.X_m(startNode), ...
                 nodes.Y_m(endNode)-nodes.Y_m(startNode)];
    globalNorm = hypot(globalVec(1),globalVec(2));

    if globalNorm<=0
        candidate = emptyCandidate();
        return;
    end

    globalUnit = globalVec/globalNorm;

    for s = 1:nStates

        e = stateEdge(s);
        dir = stateDirection(s);

        if dir==1
            stateStartNode(s) = edges.Node1(e);
            stateEndNode(s) = edges.Node2(e);
            x = edges.XGeom{e};
            y = edges.YGeom{e};
        else
            stateStartNode(s) = edges.Node2(e);
            stateEndNode(s) = edges.Node1(e);
            x = flipud(edges.XGeom{e});
            y = flipud(edges.YGeom{e});
        end

        stateStartAz(s) = localEndpointAzimuth(x,y,true);
        stateEndAz(s) = localEndpointAzimuth(x,y,false);

        % Minimum travel time.  Observed rank-1 traces are traversed at
        % reference speed; movement on junctions and artificial gaps is
        % slower by the selected multiplier.
        if edges.EdgeType(e)=="TRACE"
            edgeSpeed = W.traceSpeed;
        else
            edgeSpeed = W.offTraceSpeed;
        end

        stateTraversalCost(s) = ...
            edges.Length_m(e)/max(edgeSpeed,eps);
    end

    source = nStates+1;
    sink = nStates+2;

    S = [];
    T = [];
    C = [];

    % A BSR must leave its external start tip by travelling on an
    % observed rank-1 feature. Direct start-to-end gap solutions are not
    % admissible, even when they have lower travel time.
    startStates = find( ...
        stateStartNode==startNode & ...
        edges.EdgeType(stateEdge)=="TRACE");

    for i = 1:numel(startStates)
        s = startStates(i);
        S(end+1,1) = source; %#ok<AGROW>
        T(end+1,1) = s; %#ok<AGROW>
        C(end+1,1) = stateTraversalCost(s); %#ok<AGROW>
    end

    for a = 1:nStates

        nextStates = find(stateStartNode==stateEndNode(a));

        for ib = 1:numel(nextStates)

            b = nextStates(ib);

            if stateEdge(a)==stateEdge(b)
                continue;
            end

            turn = angularDifference(stateEndAz(a),stateStartAz(b));

            if turn>maxTurn
                continue;
            end

            transitionCost = stateTraversalCost(b);

            S(end+1,1) = a; %#ok<AGROW>
            T(end+1,1) = b; %#ok<AGROW>
            C(end+1,1) = transitionCost; %#ok<AGROW>
        end
    end

    % Likewise, the path must reach its external end tip on an observed
    % rank-1 feature. This guarantees observed support at both BSR ends.
    endStates = find( ...
        stateEndNode==endNode & ...
        edges.EdgeType(stateEdge)=="TRACE");

    for i = 1:numel(endStates)
        S(end+1,1) = endStates(i); %#ok<AGROW>
        T(end+1,1) = sink; %#ok<AGROW>
        C(end+1,1) = 0; %#ok<AGROW>
    end

    candidate = emptyCandidate();

    if isempty(startStates) || isempty(endStates) || isempty(S)
        return;
    end

    G = digraph(S,T,C,nStates+2);

    try
        [pathNodes,pathCost] = shortestpath(G,source,sink,'Method','positive');
    catch
        return;
    end

    if numel(pathNodes)<3
        return;
    end

    states = pathNodes(2:end-1);
    edgeIDs = stateEdge(states);

    if numel(unique(edgeIDs))<numel(edgeIDs)
        return;
    end

    usedTypes = edges.EdgeType(edgeIDs);
    if isempty(usedTypes) || usedTypes(1)~="TRACE" || ...
            usedTypes(end)~="TRACE" || ~any(usedTypes=="TRACE")
        return;
    end

    candidate.isValid = true;
    candidate.PathCost = pathCost;
    candidate.DirectedEdgeStates = [edgeIDs(:),stateDirection(states(:))];
end

function az = localEndpointAzimuth(x,y,isStart)

    if numel(x)<2
        az = NaN;
        return;
    end

    if isStart
        i1 = 1; i2 = 2;
    else
        i1 = numel(x)-1; i2 = numel(x);
    end

    az = mod(atan2d(x(i2)-x(i1),y(i2)-y(i1)),360);
end

function cost = computeGapCost(gap,preferredGap,W)

    if gap<=preferredGap
        cost = W.gapWeight*gap;
    else
        cost = W.gapWeight*preferredGap + ...
            W.longGapWeight*(gap-preferredGap);
    end
end

function [origin,axisVector,nodeProgress,observedExtent] = ...
    computeGlobalEndpointProgress(original,nodes)

    X = [];
    Y = [];

    for i = 1:numel(original)
        X = [X;original(i).x(:)]; %#ok<AGROW>
        Y = [Y;original(i).y(:)]; %#ok<AGROW>
    end

    origin = [mean(X,'omitnan') mean(Y,'omitnan')];
    centred = [X-origin(1),Y-origin(2)];

    covarianceMatrix = cov(centred,1);
    [V,D] = eig(covarianceMatrix);
    [~,idx] = max(diag(D));
    axisVector = V(:,idx)';
    axisVector = axisVector/max(hypot(axisVector(1),axisVector(2)),eps);

    observedProgress = centred*axisVector';
    observedExtent = max(observedProgress)-min(observedProgress);

    nodeProgress = ...
        (nodes.X_m-origin(1))*axisVector(1) + ...
        (nodes.Y_m-origin(2))*axisVector(2);
end

function candidate = emptyCandidate()

    candidate = struct( ...
        'isValid',false, ...
        'CandidateID',NaN, ...
        'CandidateRank',NaN, ...
        'Selected',false, ...
        'PathCost',Inf, ...
        'DirectedEdgeStates',zeros(0,2), ...
        'PairIndex',NaN, ...
        'WeightIndex',NaN, ...
        'StartNode',NaN, ...
        'EndNode',NaN, ...
        'Coverage',NaN, ...
        'GlobalExtentCoverage',NaN, ...
        'ReachesGlobalEndpoints',false, ...
        'NumberTraceEdges',NaN, ...
        'NumberJunctionEdges',NaN, ...
        'NumberArtificialGaps',NaN, ...
        'MaximumGap_m',NaN, ...
        'MeanGap_m',NaN, ...
        'TotalGap_m',NaN, ...
        'MeanTurn_deg',NaN, ...
        'BacktrackingFraction',NaN);
end

function candidate = evaluateDynamicCandidate( ...
    candidate,nodes,edges,totalObservedLength)

    states = candidate.DirectedEdgeStates;
    edgeIDs = states(:,1);
    directions = states(:,2);

    isTrace = edges.EdgeType(edgeIDs)=="TRACE";
    isJunction = edges.EdgeType(edgeIDs)=="JUNCTION";
    isGap = edges.EdgeType(edgeIDs)=="GAP";

    traceLength = sum(edges.Length_m(edgeIDs(isTrace)));

    candidate.Coverage = min(traceLength/max(totalObservedLength,eps),1);
    candidate.NumberTraceEdges = sum(isTrace);
    candidate.NumberJunctionEdges = sum(isJunction);
    candidate.NumberArtificialGaps = sum(isGap);

    gapLengths = edges.Length_m(edgeIDs(isGap));

    if isempty(gapLengths)
        candidate.MaximumGap_m = 0;
        candidate.MeanGap_m = 0;
        candidate.TotalGap_m = 0;
    else
        candidate.MaximumGap_m = max(gapLengths);
        candidate.MeanGap_m = mean(gapLengths);
        candidate.TotalGap_m = sum(gapLengths);
    end

    turns = [];
    backwardLength = 0;

    globalVec = [nodes.X_m(candidate.EndNode)-nodes.X_m(candidate.StartNode), ...
                 nodes.Y_m(candidate.EndNode)-nodes.Y_m(candidate.StartNode)];
    globalUnit = globalVec/max(hypot(globalVec(1),globalVec(2)),eps);

    for i = 1:numel(edgeIDs)

        [x,y,startNode,endNode] = orientedEdgeGeometry( ...
            edges,edgeIDs(i),directions(i));

        v = [nodes.X_m(endNode)-nodes.X_m(startNode), ...
             nodes.Y_m(endNode)-nodes.Y_m(startNode)];
        v = v/max(hypot(v(1),v(2)),eps);

        if dot(v,globalUnit)<0 && edges.EdgeType(edgeIDs(i))=="TRACE"
            backwardLength = backwardLength + edges.Length_m(edgeIDs(i));
        end

        if i<numel(edgeIDs)
            [xn,yn,~,~] = orientedEdgeGeometry( ...
                edges,edgeIDs(i+1),directions(i+1));
            az1 = localEndpointAzimuth(x,y,false);
            az2 = localEndpointAzimuth(xn,yn,true);
            turns(end+1,1) = angularDifference(az1,az2); %#ok<AGROW>
        end
    end

    if isempty(turns)
        candidate.MeanTurn_deg = 0;
    else
        candidate.MeanTurn_deg = mean(turns);
    end

    candidate.BacktrackingFraction = backwardLength/max(traceLength,eps);
end

function [x,y,startNode,endNode] = orientedEdgeGeometry(edges,e,direction)

    if direction==1
        x = edges.XGeom{e};
        y = edges.YGeom{e};
        startNode = edges.Node1(e);
        endNode = edges.Node2(e);
    else
        x = flipud(edges.XGeom{e});
        y = flipud(edges.YGeom{e});
        startNode = edges.Node2(e);
        endNode = edges.Node1(e);
    end
end

function d = angularDifference(a,b)

    d = abs(mod(a-b+180,360)-180);
end

function candidates = removeDuplicateDynamicCandidates(candidates)

    candidates = candidates(:);
    keys = strings(numel(candidates),1);

    for i = 1:numel(candidates)
        S = candidates(i).DirectedEdgeStates;
        token = strings(size(S,1),1);
        for j = 1:size(S,1)
            token(j) = string(S(j,1))+"_"+string(S(j,2));
        end
        keys(i) = strjoin(token,"-");
    end

    [~,ia] = unique(keys,'stable');
    candidates = candidates(ia);
end

function [selectedIndex,logText,candidates] = ...
    selectCandidateLexicographically(candidates,S)

    candidates = candidates(:);
    n = numel(candidates);

    for i = 1:n
        candidates(i).CandidateID = i;
    end

    active = (1:n)';
    logText = strings(0,1);

    % Completeness is a prerequisite whenever at least one candidate spans
    % the required fraction of the observed rupture extent. Partial paths
    % are considered only when no complete path exists at all.
    completeMask = [candidates(active).GlobalExtentCoverage]' >= ...
        S.minimumGlobalExtentForCompletion;

    if any(completeMask)
        active = active(completeMask);
        logText(end+1) = sprintf( ...
            ['1 Global completion: retained %d candidates with ', ...
             'extent >= %.3f.'], ...
            numel(active),S.minimumGlobalExtentForCompletion);
    else
        values = [candidates(active).GlobalExtentCoverage]';
        target = max(values);
        active = active(values>=target-S.globalExtentTolerance);
        logText(end+1) = sprintf( ...
            ['1 Global completion: no complete candidate; max extent ', ...
             '%.5f, tolerance %.5f, remaining %d.'], ...
            target,S.globalExtentTolerance,numel(active));
    end

    % Among complete candidates, retain those with the greatest endpoint
    % span before applying trace coverage and connection criteria.
    values = [candidates(active).GlobalExtentCoverage]';
    target = max(values);
    active = active(values>=target-S.globalExtentTolerance);
    logText(end+1) = sprintf( ...
        '2 Global extent: max %.5f, tolerance %.5f, remaining %d.', ...
        target,S.globalExtentTolerance,numel(active));

    values = [candidates(active).Coverage]';
    target = max(values);
    active = active(values>=target-S.coverageTolerance);
    logText(end+1) = sprintf( ...
        '3 Trace coverage: max %.5f, tolerance %.5f, remaining %d.', ...
        target,S.coverageTolerance,numel(active));

    values = [candidates(active).NumberArtificialGaps]';
    target = min(values);
    active = active(values<=target+S.nGapTolerance);
    logText(end+1) = sprintf( ...
        '4 Number of gaps: min %d, remaining %d.',target,numel(active));

    values = [candidates(active).MaximumGap_m]';
    target = min(values);
    active = active(values<=target+S.maxGapTolerance_m);
    logText(end+1) = sprintf( ...
        '5 Maximum gap: min %.2f m, remaining %d.',target,numel(active));

    values = [candidates(active).TotalGap_m]';
    target = min(values);
    active = active(values<=target+S.totalGapTolerance_m);
    logText(end+1) = sprintf( ...
        '6 Total gap: min %.2f m, remaining %d.',target,numel(active));

    values = [candidates(active).MeanTurn_deg]';
    target = min(values);
    active = active(values<=target+S.meanTurnTolerance);
    logText(end+1) = sprintf( ...
        '7 Mean turn: min %.3f deg, remaining %d.',target,numel(active));

    values = [candidates(active).BacktrackingFraction]';
    target = min(values);
    active = active(values<=target+S.backtrackTolerance);
    logText(end+1) = sprintf( ...
        '8 Backtracking: min %.7f, remaining %d.',target,numel(active));

    values = [candidates(active).PathCost]';
    target = min(values);
    active = active(values<=target+S.pathCostTolerance);
    logText(end+1) = sprintf( ...
        '9 Travel time: min %.6f, remaining %d.',target,numel(active));

    selectedIndex = active(1);

    sortTable = table( ...
        -[candidates.GlobalExtentCoverage]', ...
        -[candidates.Coverage]', ...
        [candidates.NumberArtificialGaps]', ...
        [candidates.MaximumGap_m]', ...
        [candidates.TotalGap_m]', ...
        [candidates.MeanTurn_deg]', ...
        [candidates.BacktrackingFraction]', ...
        [candidates.PathCost]', ...
        (1:n)', ...
        'VariableNames',{'NegGlobalExtent','NegCoverage','NGap', ...
        'MaximumGap','TotalGap','MeanTurn','Backtracking', ...
        'PathCost','OriginalIndex'});

    sortTable = sortrows(sortTable, ...
        {'NegGlobalExtent','NegCoverage','NGap','MaximumGap', ...
        'TotalGap','MeanTurn','Backtracking','PathCost'}, ...
        {'ascend','ascend','ascend','ascend','ascend','ascend', ...
         'ascend','ascend'});

    for r = 1:height(sortTable)
        candidates(sortTable.OriginalIndex(r)).CandidateRank = r;
    end

    candidates(selectedIndex).Selected = true;
    logText(end+1) = sprintf('Selected candidate ID %d.', ...
        candidates(selectedIndex).CandidateID);
end

function [edges,statesOut,nConverted] = convertSelectedJunctionsToGaps( ...
    edges,statesIn,tolerance_m)

    statesOut = statesIn;
    nConverted = 0;

    if isempty(statesIn)
        return;
    end

    for i = 1:size(statesIn,1)
        e = statesIn(i,1);

        if edges.EdgeType(e) ~= "JUNCTION"
            continue;
        end

        % A junction with a finite geometric length is an artificial
        % connection in the final BSR, although it remains a JUNCTION in
        % the internal graph used during path search.
        if edges.Length_m(e) <= tolerance_m
            continue;
        end

        newEdgeID = height(edges)+1;
        newRow = edges(e,:);
        newRow.EdgeID = newEdgeID;
        newRow.EdgeType = "GAP";

        edges = [edges;newRow]; %#ok<AGROW>
        statesOut(i,1) = newEdgeID;
        nConverted = nConverted+1;
    end
end

function [edges,statesOut,nInserted] = insertMissingAssemblyGaps( ...
    nodes,edges,statesIn,tolerance_m)

    %#ok<INUSD> nodes is retained in the signature for interface clarity.
    if isempty(statesIn) || size(statesIn,1)<2
        statesOut = statesIn;
        nInserted = 0;
        return;
    end

    statesOut = zeros(0,2);
    nInserted = 0;

    for i = 1:size(statesIn,1)
        currentState = statesIn(i,:);
        statesOut(end+1,:) = currentState; %#ok<AGROW>

        if i == size(statesIn,1)
            continue;
        end

        e1 = currentState(1);
        d1 = currentState(2);
        e2 = statesIn(i+1,1);
        d2 = statesIn(i+1,2);

        [x1,y1,~,endNode1] = orientedEdgeGeometry(edges,e1,d1);
        [x2,y2,startNode2,~] = orientedEdgeGeometry(edges,e2,d2);

        if isempty(x1) || isempty(x2)
            continue;
        end

        joinDistance_m = hypot(x1(end)-x2(1),y1(end)-y2(1));
        if joinDistance_m <= tolerance_m
            continue;
        end

        newEdgeID = height(edges)+1;
        newRow = edges(1,:);
        newRow.EdgeID = newEdgeID;
        newRow.Node1 = endNode1;
        newRow.Node2 = startNode2;
        newRow.EdgeType = "GAP";
        newRow.Length_m = joinDistance_m;
        newRow.OriginalPartIndex = NaN;
        newRow.IdS = NaN;
        newRow.StartAlong_m = NaN;
        newRow.EndAlong_m = NaN;
        newRow.XGeom = {[x1(end);x2(1)]};
        newRow.YGeom = {[y1(end);y2(1)]};

        edges = [edges;newRow]; %#ok<AGROW>
        statesOut(end+1,:) = [newEdgeID 1]; %#ok<AGROW>
        nInserted = nInserted+1;
    end
end

function BSR = assembleDynamicBSR( ...
    nodes,edges,states,lon0,lat0)

    X = [];
    Y = [];
    EdgeID = [];
    EdgeType = strings(0,1);
    OriginalPartIndex = [];
    IdS = [];

    % Logical flag for vertices belonging to connector/gap edges.
    % It must remain logical because it is later used as an index.
    IsGapVertex = false(0,1);

    % ------------------------------------------------------------------
    % Output sampling of artificial connections.
    %
    % TRACE geometries retain all their original mapped vertices. GAP
    % geometries are straight graph edges and normally contain only their
    % two endpoints. For the exported BSR vertex table, densify every GAP
    % so that distance and x/L are also sampled across artificial
    % connections at a spacing comparable with the selected TRACE
    % geometries.
    %
    % The automatic target spacing is the median positive vertex spacing
    % along selected TRACE edges. Limits prevent extremely dense or sparse
    % input traces from creating impractical GAP output tables.
    % ------------------------------------------------------------------

    minimumGapOutputStep_m = 10;
    maximumGapOutputStep_m = 100;

    traceVertexSteps = zeros(0,1);

    for iState = 1:size(states,1)
        eState = states(iState,1);
        dState = states(iState,2);

        if edges.EdgeType(eState) ~= "TRACE"
            continue;
        end

        [xt,yt,~,~] = orientedEdgeGeometry(edges,eState,dState);
        ds = hypot(diff(xt),diff(yt));
        ds = ds(isfinite(ds) & ds>0);
        traceVertexSteps = [traceVertexSteps;ds(:)]; %#ok<AGROW>
    end

    if isempty(traceVertexSteps)
        gapOutputStep_m = maximumGapOutputStep_m;
    else
        gapOutputStep_m = median(traceVertexSteps);
        gapOutputStep_m = max(minimumGapOutputStep_m, ...
            min(maximumGapOutputStep_m,gapOutputStep_m));
    end

    for i = 1:size(states,1)

        e = states(i,1);
        direction = states(i,2);
        [x,y,~,~] = orientedEdgeGeometry(edges,e,direction);

        % Densify only artificial connections. This changes neither the
        % selected path nor its geometry: all inserted points lie exactly
        % on the original straight GAP segment.
        if edges.EdgeType(e)=="GAP" && numel(x)>=2
            gapLength_m = sum(hypot(diff(x),diff(y)));
            nGapSegments = max(1,ceil(gapLength_m/gapOutputStep_m));

            if nGapSegments>1
                cumulativeGap = [0;cumsum(hypot(diff(x),diff(y)))];
                queryGap = linspace(0,cumulativeGap(end), ...
                    nGapSegments+1)';

                % Remove repeated cumulative positions defensively before
                % interpolation. GAP edges are normally two-point lines.
                [cumulativeGapUnique,ia] = unique(cumulativeGap,'stable');
                xUnique = x(ia);
                yUnique = y(ia);

                if numel(cumulativeGapUnique)>=2
                    x = interp1(cumulativeGapUnique,xUnique,queryGap,'linear');
                    y = interp1(cumulativeGapUnique,yUnique,queryGap,'linear');
                end
            end
        end

        if isempty(X)
            use = 1:numel(x);
        else
            if hypot(X(end)-x(1),Y(end)-y(1))<1e-6
                use = 2:numel(x);
            else
                use = 1:numel(x);
            end
        end

        X = [X;x(use)]; %#ok<AGROW>
        Y = [Y;y(use)]; %#ok<AGROW>
        EdgeID = [EdgeID;repmat(e,numel(use),1)]; %#ok<AGROW>
        EdgeType = [EdgeType;repmat(edges.EdgeType(e),numel(use),1)]; %#ok<AGROW>
        OriginalPartIndex = [OriginalPartIndex; ...
            repmat(edges.OriginalPartIndex(e),numel(use),1)]; %#ok<AGROW>
        IdS = [IdS;repmat(edges.IdS(e),numel(use),1)]; %#ok<AGROW>
        IsGapVertex = [IsGapVertex; ...
            repmat(edges.EdgeType(e)~="TRACE",numel(use),1)]; %#ok<AGROW>
    end

    segmentLength = hypot(diff(X),diff(Y));

    % A segment ending at a vertex that belongs to a non-TRACE edge is
    % treated as connector/gap. Force the array to logical explicitly:
    % concatenating an empty numeric array with logical values can otherwise
    % convert 0/1 flags to doubles, and MATLAB rejects 0 as an array index.
    segmentIsConnector = logical(IsGapVertex(2:end));

    segmentTraceOnly = segmentLength;
    segmentTraceOnly(segmentIsConnector) = 0;

    CumWithGaps = [0;cumsum(segmentLength)];
    CumTraceOnly = [0;cumsum(segmentTraceOnly)];

    LengthWithGaps = CumWithGaps(end);
    TraceLength = CumTraceOnly(end);

    % The exported BSR distance and x/L always include artificial gaps.
    CumUsed = CumWithGaps;
    LengthUsed = LengthWithGaps;

    [lon,lat] = local2lonlat(X,Y,lon0,lat0);

    BSR = struct();
    BSR.X = X;
    BSR.Y = Y;
    BSR.Lon = lon;
    BSR.Lat = lat;
    BSR.EdgeID = EdgeID;
    BSR.EdgeType = EdgeType;
    BSR.OriginalPartIndex = OriginalPartIndex;
    BSR.IdS = IdS;
    BSR.CumDistanceTraceOnly_m = CumTraceOnly;
    BSR.CumDistanceUsed_m = CumUsed;
    BSR.LengthWithGaps_m = LengthWithGaps;
    BSR.TraceLength_m = TraceLength;
    BSR.TotalLengthUsed_m = LengthUsed;
    BSR.XoverL_TraceOnly = CumTraceOnly/max(TraceLength,eps);
    BSR.XoverL = CumUsed/max(LengthUsed,eps);
    BSR.GapOutputStep_m = gapOutputStep_m;
end

function T = buildSelectedEdgeTable(eventID,best,edges)

    S = best.DirectedEdgeStates;
    n = size(S,1);

    EdgeID = S(:,1);
    Direction = S(:,2);

    T = table( ...
        repmat(eventID,n,1), ...
        (1:n)', ...
        EdgeID,Direction, ...
        edges.Node1(EdgeID),edges.Node2(EdgeID), ...
        edges.EdgeType(EdgeID),edges.Length_m(EdgeID), ...
        edges.OriginalPartIndex(EdgeID),edges.IdS(EdgeID), ...
        edges.StartAlong_m(EdgeID),edges.EndAlong_m(EdgeID), ...
        'VariableNames',{'IdE','PathOrder','EdgeID','Direction', ...
        'Node1','Node2','EdgeType','Length_m','OriginalPartIndex', ...
        'IdS','StartAlong_m','EndAlong_m'});
end

function T = candidatesToDynamicTable(C,eventID)

    C = C(:);
    n = numel(C);

    T = table( ...
        repmat(eventID,n,1), ...
        [C.CandidateID]',[C.CandidateRank]',[C.Selected]', ...
        [C.PairIndex]',[C.WeightIndex]',[C.StartNode]',[C.EndNode]', ...
        reshape(arrayfun(@(s)size(s.DirectedEdgeStates,1),C),[],1), ...
        [C.NumberTraceEdges]',[C.NumberJunctionEdges]', ...
        [C.GlobalExtentCoverage]',[C.ReachesGlobalEndpoints]', ...
        [C.Coverage]', ...
        [C.NumberArtificialGaps]',[C.MaximumGap_m]', ...
        [C.MeanGap_m]',[C.TotalGap_m]',[C.MeanTurn_deg]', ...
        [C.BacktrackingFraction]',[C.PathCost]', ...
        'VariableNames',{'IdE','CandidateID','CandidateRank','Selected', ...
        'PairIndex','WeightIndex','StartNode','EndNode','NEdges', ...
        'NTraceEdges','NJunctionEdges','GlobalExtentCoverage', ...
        'ReachesGlobalEndpoints','Coverage','NumberArtificialGaps','MaximumGap_m', ...
        'MeanGap_m','TotalGap_m','MeanTurn_deg','BacktrackingFraction', ...
        'PathCost'});

    T = sortrows(T,'CandidateRank','ascend');
end

function writeSelectionExplanation(filename,eventID,best,selectionLog)

    fid = fopen(filename,'w');

    if fid<0
        warning('Cannot write %s',filename);
        return;
    end

    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>

    fprintf(fid,'BACKBONE SURFACE RUPTURE - EVENT %d\n',eventID);
    fprintf(fid,'============================================================\n\n');
    fprintf(fid,'LEXICOGRAPHIC SELECTION ORDER\n');
    fprintf(fid,'1 maximum global endpoint extent\n');
    fprintf(fid,'2 maximum trace coverage\n');
    fprintf(fid,'3 minimum number of artificial connections\n');
    fprintf(fid,'4 minimum maximum connection\n');
    fprintf(fid,'5 minimum total connection length\n');
    fprintf(fid,'6 minimum mean turn\n');
    fprintf(fid,'7 minimum backtracking\n');
    fprintf(fid,'8 minimum travel time\n\n');

    fprintf(fid,'SELECTION LOG\n');
    for i = 1:numel(selectionLog)
        fprintf(fid,'%s\n',selectionLog(i));
    end

    fprintf(fid,'\nSELECTED CANDIDATE\n');
    fprintf(fid,'Candidate ID: %d\n',best.CandidateID);
    fprintf(fid,'Candidate rank: %d\n',best.CandidateRank);
    fprintf(fid,'Start node: %d\n',best.StartNode);
    fprintf(fid,'End node: %d\n',best.EndNode);
    fprintf(fid,'Global extent coverage: %.6f\n',best.GlobalExtentCoverage);
    fprintf(fid,'Reaches global endpoints: %d\n',best.ReachesGlobalEndpoints);
    fprintf(fid,'Trace coverage: %.6f\n',best.Coverage);
    fprintf(fid,'Trace edges: %d\n',best.NumberTraceEdges);
    fprintf(fid,'Artificial gaps: %d\n',best.NumberArtificialGaps);
    fprintf(fid,'Maximum gap: %.2f m\n',best.MaximumGap_m);
    fprintf(fid,'Mean gap: %.2f m\n',best.MeanGap_m);
    fprintf(fid,'Total gap: %.2f m\n',best.TotalGap_m);
    fprintf(fid,'Mean turn: %.3f deg\n',best.MeanTurn_deg);
    fprintf(fid,'Backtracking: %.7f\n',best.BacktrackingFraction);
    fprintf(fid,'Travel time: %.6f\n',best.PathCost);
end

function makeBSRDiagnosticFigure( ...
    original,nodes,edges,best,BSR,eventID,eventDir)

    fig = figure('Color','w','Position',[100 100 1250 900]);
    hold on;
    axis equal;
    box on;
    grid on;

    % Graphic hierarchy:
    %   1) observed rank-1 traces below the BSR;
    %   2) BSR trace with a wider black stroke;
    %   3) observed rank-1 traces again above the BSR;
    %   4) artificial connections and endpoint markers.
    % The second rank-1 pass keeps the observations visible, while the
    % wider BSR remains recognizable as a black outline.
    rank1Color = [0.70 0.70 0.70];
    rank1LineWidth = 1.0;
    bsrLineWidth = 2.0;
    gapLineWidth = 1.5;

    for i = 1:numel(original)
        plot(original(i).x,original(i).y,'Color',rank1Color, ...
            'LineWidth',rank1LineWidth,'HandleVisibility','off');
    end

    states = best.DirectedEdgeStates;

    hTrace = gobjects(1);
    hGap = gobjects(1);
    traceDrawn = false;
    gapDrawn = false;

    % Draw the observed parts of the BSR first.
    for i = 1:size(states,1)
        e = states(i,1);
        direction = states(i,2);
        [x,y,~,~] = orientedEdgeGeometry(edges,e,direction);

        if edges.EdgeType(e)=="TRACE"
            h = plot(x,y,'k-','LineWidth',bsrLineWidth);
            if ~traceDrawn
                hTrace = h;
                traceDrawn = true;
            end
        end
    end

    % Draw all observed rank-1 traces again above the wider BSR.
    hObserved = gobjects(1);
    for i = 1:numel(original)
        h = plot(original(i).x,original(i).y,'Color',rank1Color, ...
            'LineWidth',rank1LineWidth);
        if i==1
            hObserved = h;
        else
            set(h,'HandleVisibility','off');
        end
    end

    % Draw artificial connections last so they remain clearly visible.
    for i = 1:size(states,1)
        e = states(i,1);
        direction = states(i,2);
        [x,y,~,~] = orientedEdgeGeometry(edges,e,direction);

        if edges.EdgeType(e)=="GAP"
            h = plot(x,y,'r--','LineWidth',gapLineWidth);
            if ~gapDrawn
                hGap = h;
                gapDrawn = true;
            else
                set(h,'HandleVisibility','off');
            end
        end
    end

    hStart = scatter(nodes.X_m(best.StartNode),nodes.Y_m(best.StartNode), ...
        75,'^','filled','MarkerFaceColor',[0.2 0.7 0.2]);

    hEnd = scatter(nodes.X_m(best.EndNode),nodes.Y_m(best.EndNode), ...
        75,'v','filled','MarkerFaceColor',[0.8 0.2 0.2]);

    title(sprintf(['Event %d - Backbone Surface Rupture\n' ...
        'extent %.2f%% | trace %.2f%% | connections %d | ' ...
        'max %.0f m | confidence %s'], ...
        eventID,100*best.GlobalExtentCoverage,100*best.Coverage, ...
        best.NumberArtificialGaps, ...
        best.MaximumGap_m,best.ReconstructionConfidence));

    xlabel('Local X (m)');
    ylabel('Local Y (m)');

    handles = hObserved;
    labels = {'Observed Rank 1 traces'};

    if isgraphics(hTrace)
        handles(end+1) = hTrace;
        labels{end+1} = 'Backbone Surface Rupture - Rank 1 traces';
    end

    if isgraphics(hGap)
        handles(end+1) = hGap;
        labels{end+1} = 'Backbone Surface Rupture - Artificial connection';
    end

    handles(end+1) = hStart;
    labels{end+1} = 'Start';
    handles(end+1) = hEnd;
    labels{end+1} = 'End';

    legend(handles,labels,'Location','bestoutside');

    exportgraphics(fig,fullfile(eventDir,sprintf('BSR_IdE_%d.png',eventID)), ...
        'Resolution',220);
    savefig(fig,fullfile(eventDir,sprintf('BSR_IdE_%d.fig',eventID)));
    close(fig);
end

function makeDynamicCandidateFigures( ...
    original,nodes,edges,candidates,eventID,eventDir, ...
    lon0,lat0,nFigures)

    candidateDir = fullfile(eventDir,'CANDIDATE_FIGURES');
    if ~exist(candidateDir,'dir')
        mkdir(candidateDir);
    end

    [~,order] = sort([candidates.CandidateRank]);

    if isinf(nFigures)
        nPlot = numel(order);
    else
        nPlot = min(nFigures,numel(order));
    end

    for ii = 1:nPlot

        c = candidates(order(ii));
        B = assembleDynamicBSR(nodes,edges,c.DirectedEdgeStates, ...
            lon0,lat0);

        fig = figure('Visible','off','Color','w','Position',[100 100 1000 800]);
        hold on;
        axis equal;
        box on;
        grid on;

        for i = 1:numel(original)
            plot(original(i).x,original(i).y,'Color',[0.82 0.82 0.82], ...
                'LineWidth',0.8);
        end

        for i = 1:size(c.DirectedEdgeStates,1)
            e = c.DirectedEdgeStates(i,1);
            d = c.DirectedEdgeStates(i,2);
            [x,y,~,~] = orientedEdgeGeometry(edges,e,d);
            if edges.EdgeType(e)=="TRACE"
                plot(x,y,'k-','LineWidth',1.5);
            elseif edges.EdgeType(e)=="JUNCTION"
                continue;
            else
                plot(x,y,'r--','LineWidth',1.5);
            end
        end

        scatter(nodes.X_m(c.StartNode),nodes.Y_m(c.StartNode),90,'^','filled');
        scatter(nodes.X_m(c.EndNode),nodes.Y_m(c.EndNode),90,'v','filled');

        title(sprintf(['Event %d | candidate rank %d | ID %d\n' ...
            'coverage %.2f%% | gaps %d | max %.0f m | turn %.1f deg'], ...
            eventID,c.CandidateRank,c.CandidateID,100*c.Coverage, ...
            c.NumberArtificialGaps,c.MaximumGap_m,c.MeanTurn_deg));

        xlabel('Local X (m)');
        ylabel('Local Y (m)');

        filename = sprintf('Candidate_%03d_ID_%04d_Cov_%05.1f.png', ...
            c.CandidateRank,c.CandidateID,100*c.Coverage);

        exportgraphics(fig,fullfile(candidateDir,filename),'Resolution',180);
        close(fig);
    end
end

function makeTravelTimeScenarioFigures( ...
    original,nodes,edges,candidates,scenarios,eventID,eventDir, ...
    lon0,lat0,selectionSettings)

    scenarioDir = fullfile(eventDir,'TRAVEL_TIME_SCENARIOS');

    if ~exist(scenarioDir,'dir')
        mkdir(scenarioDir);
    end

    for iScenario = 1:numel(scenarios)

        idx = find([candidates.WeightIndex]'==iScenario);

        if isempty(idx)
            continue;
        end

        subset = candidates(idx);

        [localBestIndex,~,subset] = ...
            selectCandidateLexicographically( ...
                subset,selectionSettings);

        candidate = subset(localBestIndex);

        if isempty(candidate.DirectedEdgeStates)
            continue;
        end

        B = assembleDynamicBSR( ...
            nodes,edges,candidate.DirectedEdgeStates, ...
            lon0,lat0);

        if isempty(B.X)
            continue;
        end

        fig = figure( ...
            'Visible','off', ...
            'Color','w', ...
            'Position',[100 100 1100 850]);

        hold on;
        axis equal;
        box on;
        grid on;

        for i = 1:numel(original)
            plot(original(i).x,original(i).y, ...
                'Color',[0.80 0.80 0.80], ...
                'LineWidth',0.8);
        end

        states = candidate.DirectedEdgeStates;

        for i = 1:size(states,1)

            e = states(i,1);
            d = states(i,2);

            [x,y,~,~] = orientedEdgeGeometry(edges,e,d);

            if edges.EdgeType(e)=="TRACE"
                plot(x,y,'k-','LineWidth',1.5);
            elseif edges.EdgeType(e)=="JUNCTION"
                continue;
            else
                plot(x,y,'r--','LineWidth',1.5);
            end
        end

        scatter(nodes.X_m(candidate.StartNode), ...
            nodes.Y_m(candidate.StartNode),100,'^','filled');

        scatter(nodes.X_m(candidate.EndNode), ...
            nodes.Y_m(candidate.EndNode),100,'v','filled');

        multiplier = scenarios(iScenario).offTraceMultiplier;

        title({ ...
            sprintf('Event %d - Minimum travel-time BSR',eventID), ...
            sprintf(['off-trace time multiplier x%d | coverage %.1f%% | ' ...
                     'gaps %d | max gap %.0f m | travel time %.0f'], ...
                multiplier,100*candidate.Coverage, ...
                candidate.NumberArtificialGaps, ...
                candidate.MaximumGap_m,candidate.PathCost)});

        xlabel('Local X (m)');
        ylabel('Local Y (m)');

        filename = sprintf( ...
            'BSR_IdE_%d_TravelTime_x%d.png', ...
            eventID,multiplier);

        exportgraphics(fig,fullfile(scenarioDir,filename), ...
            'Resolution',200);

        close(fig);
    end
end

function [confidence,geometryFlag] = ...
    assignReconstructionConfidence( ...
    globalExtentCoverage,traceCoverage,fractionOnRank1, ...
    maximumConnection_m,backtrackingFraction)

    if globalExtentCoverage<0.90 || backtrackingFraction>0.05
        confidence = "D";

    elseif maximumConnection_m>10000 || ...
            globalExtentCoverage<0.95 || fractionOnRank1<0.60
        confidence = "C";

    elseif maximumConnection_m>5000 || fractionOnRank1<0.80
        confidence = "B";

    else
        confidence = "A";
    end

    if traceCoverage<0.50
        geometryFlag = "LOW_TRACE_COVERAGE";
    else
        geometryFlag = "OK";
    end
end

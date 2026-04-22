function TE_bits = conditional_TE_delay(X, Y, Z, L, delay, deltaZ, phase)
% phase: 0/1 或者直接传 spec 结构体（推荐）

    % defaults
    if nargin < 5 || isempty(delay),    delay    = 1; end
    if nargin < 6 || isempty(deltaZ), deltaZ = 1; end
    if nargin < 7 || isempty(phase),  phase  = 0; end

    % phase -> spec（保持你现在的调用方式不变）
    if isstruct(phase)
        specX = phase; specY = phase; specZ = phase;
    else
        if phase == 1
            specX = struct('mode','phase');
            specY = struct('mode','phase');
            specZ = struct('mode','phase');
        else
            specX = struct('mode','linear');
            specY = struct('mode','linear');
            specZ = struct('mode','linear');
        end
    end

    % sanity
    if delay < 1, error('tau must be >=1'); end
    if deltaZ < 0, error('deltaZ must be >=0'); end

    N = size(Y,1);
    if size(X,1) ~= N, error('X and Y length mismatch'); end
    if ~isempty(Z) && size(Z,1) ~= N, error('Z and Y length mismatch'); end

    t_idx = (L + max([1, delay, deltaZ])) : N;
    if isempty(t_idx)
        TE_bits = NaN;
        warning('有效样本不足（L=%d, tau=%d, N=%d）', L, delay, N);
        return;
    end

    % --- build embeddings on RAW space ---
%     Y_now_raw = Y(t_idx + delay, :);
%     X_past_raw   = build_embedding(X, t_idx, L, 0);
%     Y_past_raw   = build_embedding(Y, t_idx, L, 0);
%     X_past_raw = build_embedding(X, t_idx, L, 1);
%     Y_past_raw = build_embedding(Y, t_idx, L, 1);
% 
    Y_now_raw  = Y(t_idx, :);
    X_past_raw = build_embedding(X, t_idx, L, delay);
    Y_past_raw = build_embedding(Y, t_idx, L, 1);

    if isempty(Z)
        Z_past_raw = [];
    else
        Z_past_raw = build_embedding(Z, t_idx, L, deltaZ);
    end

    % --- prep AFTER embedding (key change) ---
    Y_future = prep_vars(Y_now_raw, specY);
    X_past   = prep_vars(X_past_raw,   specX);
    Y_past   = prep_vars(Y_past_raw,   specY);

    if isempty(Z)
        condition = Y_past;
    else
        Z_past   = prep_vars(Z_past_raw, specZ);
        condition = [Y_past, Z_past];
    end

    % TE via Gaussian-copula CMI (internal copnorm)
    TE_bits = gccmi_ccc(Y_future, X_past, condition);
end


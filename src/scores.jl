"""
    fuzzy_scores(predicted, real; neighborhood, ignored=())

Return a fuzzy score for each cell. Compares each cell of the two arrays.
It finds the shortest distance within the neighborhood where the
value in the `real` array is present in `predict` array.

The longer the distance, the lower score for this cell. Maximum
is 1 and if not present its 0. The weight function is linear.

- `predict`: Array of predicted data
- `real`: Array of real target data

# Keywords

- `neighborhood`: any `Neighborhood` oject, `Window(1)` by default.
- `banlist`: Values in the arrays which will not be evaluated.

Adapted from: https://github.com/johanlahti/urban-lu-model/blob/master/script/Utils.py
"""
function fuzzy_scores(predicted, real; neighborhood=Window{1}(), ignored=0)
    broadcast(real, predicted, CartesianIndices(real)) do r, p, I
        r in ignored && return 0.0

        weight = 0.0
        shortestdist = Inf

        # Get neighborhood values from the `real` values array
        real_hood = _setbuffer(neighborhood, _getwindow(real, Tuple(I)))
        # Search the neighborhood
        for (hood_val, dist) in zip(neighbors(real_hood), distances(real_hood))
            if p === hood_val && dist < shortestdist
                # Store the shortest distance at which we found the value
                shortestdist = dist
                weight = 1 - dist / (length(hood) + 1)
            end
        end
        return weight
    end
end

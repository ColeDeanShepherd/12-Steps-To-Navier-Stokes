#include "GraphMetrics.h"

Vector2d graphPointToPxPoint(const GraphMetrics& graphMetrics, const Vector2d& graphPoint)
{
    const auto xPercentFromLeft = (graphPoint.x - graphMetrics.minX) / (graphMetrics.maxX - graphMetrics.minX);
    const auto yPercentFromBottom = (graphPoint.y - graphMetrics.minY) / (graphMetrics.maxY - graphMetrics.minY);

    return Vector2d(
        graphMetrics.pos.x + (xPercentFromLeft * graphMetrics.width),
        (graphMetrics.pos.y + graphMetrics.height) - (yPercentFromBottom * graphMetrics.height)
    );
}
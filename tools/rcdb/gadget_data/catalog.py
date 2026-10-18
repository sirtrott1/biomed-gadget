"""Paged administrative indexes over immutable native record references."""
from bisect import bisect_left
from collections.abc import Sequence
from copy import deepcopy


class CoordinateSpec:
    def __init__(self, parts):
        self.parts = parts

    def __len__(self):
        return sum(map(len, self.parts))

    def __iter__(self):
        for part in self.parts:
            yield from part


class WriteIndex:
    """Keep only pending references while publishing immutable index pages."""
    def __init__(self, writer, kind, name):
        self.writer, self.kind, self.name = writer, kind, name
        self.pending, self.pages, self.weight = [], [], 0
        self.count, self.extent, self.area, self.occurrences = 0, 0, 0, 0
        self.last = None

    def __len__(self):
        return self.count

    def __getitem__(self, index):
        if index != -1 or self.last is None:
            raise IndexError(index)
        return self.last

    def append(self, reference):
        from rexgraph.value_codec import pack_value
        weight = len(pack_value(reference))
        if self.pending and (len(self.pending) == 128 or self.weight + weight >
                             min(self.writer.limits.record_bytes, self.writer.limits.decode_bytes) // 2):
            self.flush()
        self.pending.append(reference)
        self.weight += weight
        self.count += 1
        self.last = reference
        if self.kind == "axis":
            self.extent = reference["offset"] + reference["count"]
        elif self.kind == "field":
            self.area += reference["shape"][0] * reference["shape"][1]
            self.occurrences += reference.get("entry_count", 0)

    def flush(self):
        self.pages.extend(self.writer._pages(self.kind, self.name, self.pending, start_index=len(self.pages)))
        self.pending, self.weight = [], 0

    def finish(self):
        self.flush()
        return self.pages


def intersects(selected, lower, upper):
    if selected is None:
        return True
    position = bisect_left(selected, lower)
    return position < len(selected) and selected[position] < upper


class Entries(Sequence):
    def __init__(self, store, pages):
        self.store, self.pages = store, tuple(pages)
        self.count = sum(page["count"] if "page" in page else 1 for page in self.pages)

    def __len__(self):
        return self.count

    def _entries(self, page):
        return self.store._page(page["page"]) if "page" in page else (page,)

    def __iter__(self):
        for page in self.pages:
            for reference in self._entries(page):
                yield deepcopy(reference)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(self)[index]
        if index < 0:
            index += self.count
        if not 0 <= index < self.count:
            raise IndexError(index)
        for page in self.pages:
            count = page["count"] if "page" in page else 1
            if index < count:
                return deepcopy(self._entries(page)[index])
            index -= count
        raise IndexError(index)

    @property
    def extent(self):
        return max((page["extent"][1] if "page" in page else page["offset"] + page["count"]
                    for page in self.pages), default=0)

    def axis(self, selected):
        selected = sorted(selected)
        for page in self.pages:
            lower, upper = page["extent"] if "page" in page else (page["offset"], page["offset"] + page["count"])
            if intersects(selected, lower, upper):
                for block in self._entries(page):
                    if intersects(selected, block["offset"], block["offset"] + block["count"]):
                        yield deepcopy(block)

    def field(self, rows, columns):
        rows = None if rows is None else sorted(rows)
        columns = None if columns is None else sorted(columns)
        for page in self.pages:
            if "page" in page:
                lower, upper = page["bounds"]
            else:
                lower = page["offset"]
                upper = tuple(a + b for a, b in zip(lower, page["shape"], strict=True))
            if not intersects(rows, lower[0], upper[0]) or not intersects(columns, lower[1], upper[1]):
                continue
            for block in self._entries(page):
                lower = block["offset"]
                upper = tuple(a + b for a, b in zip(lower, block["shape"], strict=True))
                if intersects(rows, lower[0], upper[0]) and intersects(columns, lower[1], upper[1]):
                    yield deepcopy(block)

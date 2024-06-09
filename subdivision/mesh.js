class HalfEdge {
    v
    next
    twin
    constructor(v, next = null, twin = null, face = null) {
        this.v = v;
        this.next = next;
        this.twin = twin;
    }

    sideCount() {
        let ans = 1
        for(let he=this.next; he!=this; he=he.next) ans += 1
        return ans
    }
}
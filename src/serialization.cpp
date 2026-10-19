#include "detail.hpp"
#include <cstring>

namespace fxp_physics {
using namespace detail;
namespace {
constexpr std::uint32_t format_version = 1;
constexpr std::uint64_t math_version = UINT64_C(0xf18e16b2f6c404da);
constexpr std::size_t max_save_bytes = 256u * 1024u * 1024u;
std::uint64_t hash_bytes(const Bytes &bytes, std::size_t count) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (std::size_t i = 0; i < count; ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}
struct Writer {
    Bytes bytes;
    void u8(std::uint8_t v) { bytes.push_back(v); }
    void u32(std::uint32_t v) {
        for (unsigned i = 0; i < 4; ++i)
            u8(static_cast<std::uint8_t>(v >> (i * 8)));
    }
    void u64(std::uint64_t v) {
        for (unsigned i = 0; i < 8; ++i)
            u8(static_cast<std::uint8_t>(v >> (i * 8)));
    }
    void real(Real v) { u64(static_cast<std::uint64_t>(v.raw_value())); }
    void boolean(bool v) { u8(v ? 1 : 0); }
    void vec(Vec3 v) {
        real(v.x);
        real(v.y);
        real(v.z);
    }
    void pose(Pose p) {
        vec(p.position);
        const auto q = p.rotation.get_internal_vector();
        real(q.x);
        real(q.y);
        real(q.z);
        real(q.w);
    }
    void inertia(Inertia m) {
        real(m.xx);
        real(m.yy);
        real(m.zz);
        real(m.xy);
        real(m.xz);
        real(m.yz);
    }
    void material(Material m) {
        real(m.static_friction);
        real(m.dynamic_friction);
        real(m.restitution);
    }
    void count(std::size_t n) {
        require(n <= 1000000, "save element limit exceeded");
        u32(static_cast<std::uint32_t>(n));
    }
};
struct Reader {
    const Bytes &bytes;
    std::size_t cursor = 0, end;
    explicit Reader(const Bytes &b) : bytes(b), end(b.size() - 8) {}
    std::uint8_t u8() {
        require(cursor < end, "truncated save");
        return bytes[cursor++];
    }
    std::uint32_t u32() {
        std::uint32_t v = 0;
        for (unsigned i = 0; i < 4; ++i)
            v |= std::uint32_t(u8()) << (i * 8);
        return v;
    }
    std::uint64_t u64() {
        std::uint64_t v = 0;
        for (unsigned i = 0; i < 8; ++i)
            v |= std::uint64_t(u8()) << (i * 8);
        return v;
    }
    Real real() {
        const auto v = u64();
        return Real::from_raw(v <= std::uint64_t(INT64_MAX) ? static_cast<std::int64_t>(v)
                                                            : -1 - static_cast<std::int64_t>(~v));
    }
    bool boolean() {
        const auto v = u8();
        require(v <= 1, "invalid boolean");
        return v != 0;
    }
    Vec3 vec() {
        const Real x = real(), y = real(), z = real();
        return {x, y, z};
    }
    Pose pose() {
        Pose p;
        p.position = vec();
        const Real x = real(), y = real(), z = real(), w = real();
        p.rotation = Quat(fxp::vec4<Real>{x, y, z, w});
        return p;
    }
    Inertia inertia() {
        Inertia m;
        m.xx = real();
        m.yy = real();
        m.zz = real();
        m.xy = real();
        m.xz = real();
        m.yz = real();
        return m;
    }
    Material material() {
        Material m;
        m.static_friction = real();
        m.dynamic_friction = real();
        m.restitution = real();
        return m;
    }
    std::uint32_t count(std::size_t minimum_bytes = 1) {
        const auto n = u32();
        if (n > 1000000 || n > (end - cursor) / minimum_bytes)
            throw std::invalid_argument("invalid save count at byte " + std::to_string(cursor) + ", count " +
                                        std::to_string(n) + ", minimum " + std::to_string(minimum_bytes) +
                                        ", remaining " + std::to_string(end - cursor));
        return n;
    }
};
void write_config(Writer &w, const Config &c) {
    w.real(c.timestep);
    w.u32(c.substeps);
    w.u32(c.velocity_iterations);
    w.vec(c.gravity);
    w.real(c.contact_slop);
    w.real(c.contact_margin);
    w.real(c.stabilization);
    w.real(c.max_correction_speed);
    w.real(c.restitution_threshold);
    w.real(c.sleep_linear_speed);
    w.real(c.sleep_angular_speed);
    w.u32(c.sleep_ticks);
    w.real(c.max_linear_speed);
    w.real(c.max_angular_speed);
}
Config read_config(Reader &r) {
    Config c;
    c.timestep = r.real();
    c.substeps = r.u32();
    c.velocity_iterations = r.u32();
    c.gravity = r.vec();
    c.contact_slop = r.real();
    c.contact_margin = r.real();
    c.stabilization = r.real();
    c.max_correction_speed = r.real();
    c.restitution_threshold = r.real();
    c.sleep_linear_speed = r.real();
    c.sleep_angular_speed = r.real();
    c.sleep_ticks = r.u32();
    c.max_linear_speed = r.real();
    c.max_angular_speed = r.real();
    validate_config(c);
    return c;
}
void write_shape(Writer &w, const Shape &shape) {
    const auto &s = shape.data();
    w.u8(static_cast<std::uint8_t>(s.type));
    switch (s.type) {
    case ShapeType::box:
        w.vec(s.half_extents);
        break;
    case ShapeType::sphere:
        w.real(s.radius);
        break;
    case ShapeType::capsule:
        w.real(s.radius);
        w.real(s.half_segment);
        break;
    case ShapeType::convex:
        w.count(s.vertices.size());
        for (auto v : s.vertices)
            w.vec(v);
        w.count(s.faces.size());
        for (const auto &f : s.faces) {
            w.count(f.size());
            for (auto i : f)
                w.u32(i);
        }
        break;
    case ShapeType::mesh:
        w.count(s.vertices.size());
        for (auto v : s.vertices)
            w.vec(v);
        w.count(s.mesh_triangles.size());
        for (auto t : s.mesh_triangles) {
            w.u32(t.a);
            w.u32(t.b);
            w.u32(t.c);
        }
        break;
    case ShapeType::heightfield:
        w.u32(s.width);
        w.u32(s.depth);
        w.real(s.spacing);
        w.count(s.heights.size());
        for (auto h : s.heights)
            w.real(h);
        break;
    }
}
ShapePtr read_shape(Reader &r) {
    const auto type = static_cast<ShapeType>(r.u8());
    switch (type) {
    case ShapeType::box:
        return Shape::box(r.vec());
    case ShapeType::sphere:
        return Shape::sphere(r.real());
    case ShapeType::capsule: {
        const Real radius = r.real(), half = r.real();
        return Shape::capsule(radius, half);
    }
    case ShapeType::convex:
    case ShapeType::mesh: {
        const auto count = r.count(24);
        std::vector<Vec3> vertices;
        vertices.reserve(count);
        for (std::uint32_t i = 0; i < count; ++i)
            vertices.push_back(r.vec());
        const auto face_count = r.count(4);
        if (type == ShapeType::convex) {
            require(face_count <= 256, "too many hull faces");
            std::vector<std::vector<std::uint32_t>> faces;
            for (std::uint32_t i = 0; i < face_count; ++i) {
                const auto n = r.count(4);
                require(n <= 128, "too many face vertices");
                std::vector<std::uint32_t> f;
                for (std::uint32_t j = 0; j < n; ++j)
                    f.push_back(r.u32());
                faces.push_back(std::move(f));
            }
            return Shape::convex(std::move(vertices), std::move(faces));
        }
        std::vector<Triangle> triangles;
        triangles.reserve(face_count);
        for (std::uint32_t i = 0; i < face_count; ++i) {
            Triangle t;
            t.a = r.u32();
            t.b = r.u32();
            t.c = r.u32();
            triangles.push_back(t);
        }
        return Shape::mesh(std::move(vertices), std::move(triangles));
    }
    case ShapeType::heightfield: {
        const auto width = r.u32(), depth = r.u32();
        const Real spacing = r.real();
        const auto n = r.count(8);
        std::vector<Real> heights;
        heights.reserve(n);
        for (std::uint32_t i = 0; i < n; ++i)
            heights.push_back(r.real());
        return Shape::heightfield(width, depth, spacing, std::move(heights));
    }
    }
    throw std::invalid_argument("unknown shape tag");
}
void validate_saved_pose(Pose p, Real position_limit = 8192_r) {
    (void)validated_pose(p, position_limit);
    require(fxp::abs(p.rotation.dot(p.rotation) - 1_r) < 0.00001_r, "save rotation is not unit length");
}
void validate_unit(Vec3 v) {
    require(fxp::abs(length(v) - 1_r) < 0.00001_r, "save normal is not unit length");
}
} // namespace

Bytes World::save() const {
    const auto &s = *impl_;
    if (s.stepping)
        throw std::logic_error("save requires a completed tick boundary");
    Writer w;
    for (const auto c : std::array<std::uint8_t, 8>{'F', 'X', 'P', 'P', 'H', 'Y', 'S', 0})
        w.u8(c);
    w.u32(format_version);
    w.u64(math_version);
    write_config(w, s.config);
    w.u64(s.tick);
    w.u64(s.rng.state());
    w.u64(s.next_body);
    w.u64(s.next_collider);
    w.u64(s.next_spring);
    for (auto layer : s.layers)
        w.u32(layer);
    w.count(s.bodies.size());
    for (const auto &b : s.bodies) {
        const auto &d = b.desc;
        w.u64(b.id);
        w.u8(static_cast<std::uint8_t>(d.type));
        w.pose(d.pose);
        w.vec(d.linear_velocity);
        w.vec(d.angular_velocity);
        w.real(d.mass);
        w.vec(d.center_of_mass);
        w.inertia(d.inertia);
        w.real(d.linear_damping);
        w.real(d.angular_damping);
        w.real(d.gravity_scale);
        w.u8(d.rotation_locks);
        w.boolean(d.allow_sleep);
        w.u64(d.user_data);
        w.pose(b.previous);
        w.pose(b.target);
        w.boolean(b.has_target);
        w.boolean(b.sleeping);
        w.u32(b.quiet_ticks);
        w.vec(b.force);
        w.vec(b.torque);
        w.real(b.inverse_mass);
        w.inertia(b.inverse_inertia);
    }
    // Serialize assets by value: sharing/allocation choices cannot change bytes.
    w.count(s.colliders.size());
    for (const auto &c : s.colliders) {
        w.u64(c.id);
        w.u64(c.desc.body);
        write_shape(w, *c.desc.shape);
        w.pose(c.desc.local_pose);
        w.u64(c.desc.parent);
        w.material(c.desc.material);
        w.u8(c.desc.layer);
        w.boolean(c.desc.sensor);
        w.u64(c.desc.user_data);
        require(w.bytes.size() <= max_save_bytes, "save exceeds 256 MiB limit");
    }
    w.count(s.springs.size());
    for (const auto &spring : s.springs) {
        w.u64(spring.id);
        w.u64(spring.desc.body_a);
        w.u64(spring.desc.body_b);
        w.vec(spring.desc.anchor_a);
        w.vec(spring.desc.anchor_b);
        w.real(spring.desc.rest_length);
        w.real(spring.desc.frequency);
        w.real(spring.desc.damping_ratio);
        w.real(spring.desc.max_force);
        w.boolean(spring.desc.collide_connected);
        w.real(spring.impulse);
    }
    w.count(s.contacts.size());
    for (const auto &c : s.contacts) {
        w.u64(c.a);
        w.u64(c.b);
        w.u32(c.triangle_a);
        w.u32(c.triangle_b);
        w.vec(c.local_a);
        w.vec(c.local_b);
        w.vec(c.normal);
        w.vec(c.tangent_impulse);
        w.real(c.normal_impulse);
        w.real(c.separation);
        w.real(c.bounce);
    }
    w.count(s.pairs.size());
    for (const auto &p : s.pairs) {
        w.u64(p.a);
        w.u64(p.b);
        w.u64(p.body_a);
        w.u64(p.body_b);
        w.boolean(p.sensor);
        w.vec(p.point);
        w.vec(p.normal);
    }
    w.count(s.events.size());
    for (const auto &e : s.events) {
        w.u64(e.tick);
        w.u64(e.collider_a);
        w.u64(e.collider_b);
        w.u64(e.body_a);
        w.u64(e.body_b);
        w.u8(static_cast<std::uint8_t>(e.type));
        w.boolean(e.sensor);
        w.vec(e.point);
        w.vec(e.normal);
    }
    w.u32(s.stats.broadphase_pairs);
    w.u32(s.stats.contacts);
    w.u64(hash_bytes(w.bytes, w.bytes.size()));
    require(w.bytes.size() <= max_save_bytes, "save exceeds 256 MiB limit");
    return std::move(w.bytes);
}
std::uint64_t World::checksum() const {
    const auto bytes = save();
    return hash_bytes(bytes, bytes.size() - 8);
}

void World::load(const Bytes &bytes) {
    impl_->mutable_access();
    if (impl_->stepping)
        throw std::logic_error("load requires a completed tick boundary");
    require(bytes.size() >= 32 && bytes.size() <= max_save_bytes, "invalid save size");
    std::uint64_t stored = 0;
    for (unsigned i = 0; i < 8; ++i)
        stored |= std::uint64_t(bytes[bytes.size() - 8 + i]) << (i * 8);
    require(stored == hash_bytes(bytes, bytes.size() - 8), "save checksum mismatch");
    Reader r(bytes);
    for (const auto c : std::array<std::uint8_t, 8>{'F', 'X', 'P', 'P', 'H', 'Y', 'S', 0})
        require(r.u8() == c, "invalid save magic");
    require(r.u32() == format_version && r.u64() == math_version, "incompatible simulation/save version");
    auto next = std::make_unique<Impl>(read_config(r), 0);
    auto &s = *next;
    s.tick = r.u64();
    s.rng.restore(r.u64());
    s.next_body = r.u64();
    s.next_collider = r.u64();
    s.next_spring = r.u64();
    require(s.next_body > 0 && s.next_collider > 0 && s.next_spring > 0, "invalid next ID");
    for (auto &layer : s.layers)
        layer = r.u32();
    for (unsigned a = 0; a < 32; ++a)
        for (unsigned b = 0; b < 32; ++b)
            require(((s.layers[a] >> b) & 1u) == ((s.layers[b] >> a) & 1u), "asymmetric layer matrix");
    const auto body_count = r.count(400);
    s.bodies.reserve(body_count);
    for (std::uint32_t i = 0; i < body_count; ++i) {
        Impl::Body b;
        auto &d = b.desc;
        b.id = r.u64();
        d.type = static_cast<BodyType>(r.u8());
        d.pose = r.pose();
        d.linear_velocity = r.vec();
        d.angular_velocity = r.vec();
        d.mass = r.real();
        d.center_of_mass = r.vec();
        d.inertia = r.inertia();
        d.linear_damping = r.real();
        d.angular_damping = r.real();
        d.gravity_scale = r.real();
        d.rotation_locks = r.u8();
        d.allow_sleep = r.boolean();
        d.user_data = r.u64();
        b.previous = r.pose();
        b.target = r.pose();
        b.has_target = r.boolean();
        b.sleeping = r.boolean();
        b.quiet_ticks = r.u32();
        b.force = r.vec();
        b.torque = r.vec();
        b.inverse_mass = r.real();
        b.inverse_inertia = r.inertia();
        require(b.id > 0 && b.id < s.next_body && (s.bodies.empty() || b.id > s.bodies.back().id),
                "invalid body ID order");
        require(static_cast<unsigned>(d.type) <= 2 && d.mass >= 0.01_r && d.mass <= 10000,
                "invalid saved body");
        validate_saved_pose(d.pose);
        validate_saved_pose(b.previous);
        validate_saved_pose(b.target);
        validate_vec(d.center_of_mass, 128, "invalid center of mass");
        validate_vec(d.linear_velocity, 128, "invalid velocity");
        validate_vec(d.angular_velocity, 32, "invalid angular velocity");
        validate_vec(b.force, 1000000, "invalid force");
        validate_vec(b.torque, 1000000, "invalid torque");
        require(d.linear_damping >= 0 && d.linear_damping <= 100 && d.angular_damping >= 0 &&
                    d.angular_damping <= 100 && fxp::abs(d.gravity_scale) <= 10 && d.rotation_locks <= 7,
                "invalid body settings");
        require(b.quiet_ticks <= s.config.sleep_ticks &&
                    (!b.sleeping ||
                     (d.type == BodyType::dynamic && zero(d.linear_velocity) && zero(d.angular_velocity))) &&
                    (!b.has_target || d.type == BodyType::kinematic),
                "invalid body lifecycle state");
        require(d.type != BodyType::static_body || (zero(d.linear_velocity) && zero(d.angular_velocity)),
                "moving static body");
        const Real expected_mass = d.type == BodyType::dynamic ? 1_r / d.mass : 0_r;
        require(b.inverse_mass == expected_mass, "invalid inverse mass");
        if (d.inertia.xx != 0 || d.inertia.yy != 0 || d.inertia.zz != 0 || d.inertia.xy != 0 ||
            d.inertia.xz != 0 || d.inertia.yz != 0)
            (void)invert(d.inertia);
        s.bodies.push_back(b);
    }
    const auto collider_count = r.count(80);
    s.colliders.reserve(collider_count);
    for (std::uint32_t i = 0; i < collider_count; ++i) {
        Impl::Collider c;
        c.id = r.u64();
        c.desc.body = r.u64();
        c.desc.shape = read_shape(r);
        c.desc.local_pose = r.pose();
        c.desc.parent = r.u64();
        c.desc.material = r.material();
        c.desc.layer = r.u8();
        c.desc.sensor = r.boolean();
        c.desc.user_data = r.u64();
        require(c.id > 0 && c.id < s.next_collider && (s.colliders.empty() || c.id > s.colliders.back().id),
                "invalid collider ID order");
        const auto &body = s.get(c.desc.body);
        validate_saved_pose(c.desc.local_pose, 128);
        validate_material(c.desc.material);
        require(c.desc.layer < 32, "invalid layer");
        const auto type = c.desc.shape->type();
        require((type != ShapeType::mesh && type != ShapeType::heightfield) ||
                    body.desc.type == BodyType::static_body,
                "moving concave shape");
        s.colliders.push_back(std::move(c));
    }
    for (const auto &c : s.colliders) {
        if (c.desc.parent)
            require(s.colliders[s.collider_index(c.desc.parent)].desc.body == c.desc.body,
                    "invalid parent body");
        validate_vec(s.local_pose(c).position, 128, "invalid compound extent");
    }
    // Inverse inertia is simulation state, but must agree exactly with its inputs.
    for (auto &b : s.bodies) {
        const Inertia expected = b.inverse_inertia;
        s.recalculate_inertia(b.id);
        const auto actual = b.inverse_inertia;
        require(expected.xx == actual.xx && expected.yy == actual.yy && expected.zz == actual.zz &&
                    expected.xy == actual.xy && expected.xz == actual.xz && expected.yz == actual.yz,
                "invalid inverse inertia");
    }
    const auto spring_count = r.count(100);
    s.springs.reserve(spring_count);
    for (std::uint32_t i = 0; i < spring_count; ++i) {
        Impl::Spring spring;
        auto &d = spring.desc;
        spring.id = r.u64();
        d.body_a = r.u64();
        d.body_b = r.u64();
        d.anchor_a = r.vec();
        d.anchor_b = r.vec();
        d.rest_length = r.real();
        d.frequency = r.real();
        d.damping_ratio = r.real();
        d.max_force = r.real();
        d.collide_connected = r.boolean();
        spring.impulse = r.real();
        require(spring.id > 0 && spring.id < s.next_spring &&
                    (s.springs.empty() || spring.id > s.springs.back().id),
                "invalid spring order");
        s.get(d.body_a);
        if (d.body_b)
            s.get(d.body_b);
        require(d.body_a != d.body_b, "self spring");
        validate_vec(d.anchor_a, 128, "invalid spring anchor");
        validate_vec(d.anchor_b, d.body_b ? 128_r : 8192_r, "invalid spring anchor");
        require(d.rest_length >= 0 && d.rest_length <= 256 && d.frequency > 0 && d.frequency <= 30 &&
                    d.damping_ratio >= 0 && d.damping_ratio <= 4 && d.max_force > 0 && d.max_force <= 1000000,
                "invalid spring settings");
        require(fxp::abs(spring.impulse) <= d.max_force * s.config.timestep, "invalid spring impulse");
        s.springs.push_back(spring);
    }
    const auto contact_count = r.count(144);
    s.contacts.reserve(contact_count);
    const auto validate_triangle = [](const Shape &shape, std::uint32_t triangle) {
        const auto &data = shape.data();
        if (data.type == ShapeType::mesh)
            require(triangle < data.mesh_triangles.size(), "invalid contact triangle");
        else if (data.type == ShapeType::heightfield)
            require(triangle < (data.width - 1) * (data.depth - 1) * 2, "invalid terrain triangle");
        else
            require(triangle == 0xffffffffu, "triangle on primitive");
    };
    for (std::uint32_t i = 0; i < contact_count; ++i) {
        Impl::Contact c;
        c.a = r.u64();
        c.b = r.u64();
        c.triangle_a = r.u32();
        c.triangle_b = r.u32();
        c.local_a = r.vec();
        c.local_b = r.vec();
        c.normal = r.vec();
        c.tangent_impulse = r.vec();
        c.normal_impulse = r.real();
        c.separation = r.real();
        c.bounce = r.real();
        require(c.a < c.b &&
                    (s.contacts.empty() ||
                     std::make_pair(s.contacts.back().a, s.contacts.back().b) <= std::make_pair(c.a, c.b)),
                "invalid contact order");
        const auto &a = s.colliders[s.collider_index(c.a)];
        const auto &b = s.colliders[s.collider_index(c.b)];
        require(a.desc.body != b.desc.body && !a.desc.sensor && !b.desc.sensor, "invalid contact owners");
        validate_triangle(*a.desc.shape, c.triangle_a);
        validate_triangle(*b.desc.shape, c.triangle_b);
        validate_vec(c.local_a, 16384, "invalid contact point");
        validate_vec(c.local_b, 16384, "invalid contact point");
        validate_unit(c.normal);
        require(c.normal_impulse >= 0 && c.normal_impulse <= 100000000 && fxp::abs(c.separation) <= 16384 &&
                    c.bounce >= 0 && c.bounce <= 256,
                "invalid contact values");
        validate_vec(c.tangent_impulse, 400000000, "invalid friction impulse");
        s.contacts.push_back(c);
    }
    const auto pair_count = r.count(81);
    s.pairs.reserve(pair_count);
    for (std::uint32_t i = 0; i < pair_count; ++i) {
        Impl::Pair p;
        p.a = r.u64();
        p.b = r.u64();
        p.body_a = r.u64();
        p.body_b = r.u64();
        p.sensor = r.boolean();
        p.point = r.vec();
        p.normal = r.vec();
        require(p.a > 0 && p.a < p.b && p.b < s.next_collider && p.body_a > 0 && p.body_b > 0 &&
                    p.body_a < s.next_body && p.body_b < s.next_body,
                "invalid pair IDs");
        require(s.pairs.empty() ||
                    std::make_pair(s.pairs.back().a, s.pairs.back().b) < std::make_pair(p.a, p.b),
                "invalid pair order");
        validate_vec(p.point, 16384, "invalid pair point");
        validate_unit(p.normal);
        s.pairs.push_back(p);
    }
    const auto event_count = r.count(90);
    s.events.reserve(event_count);
    for (std::uint32_t i = 0; i < event_count; ++i) {
        CollisionEvent e;
        e.tick = r.u64();
        e.collider_a = r.u64();
        e.collider_b = r.u64();
        e.body_a = r.u64();
        e.body_b = r.u64();
        e.type = static_cast<EventType>(r.u8());
        e.sensor = r.boolean();
        e.point = r.vec();
        e.normal = r.vec();
        require(e.tick == s.tick && static_cast<unsigned>(e.type) <= 2 && e.collider_a > 0 &&
                    e.collider_a < e.collider_b && e.collider_b < s.next_collider && e.body_a > 0 &&
                    e.body_b > 0 && e.body_a < s.next_body && e.body_b < s.next_body,
                "invalid event");
        validate_vec(e.point, 16384, "invalid event point");
        validate_unit(e.normal);
        s.events.push_back(e);
    }
    s.stats.broadphase_pairs = r.u32();
    s.stats.contacts = r.u32();
    require(r.cursor == r.end, "trailing save bytes");
    // Every check precedes the swap: malformed data cannot partially replace a world.
    impl_.swap(next);
}
} // namespace fxp_physics

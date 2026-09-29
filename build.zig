const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const bitdp = b.dependency("bitdp", .{ .target = target, .optimize = optimize }).module("bitdp");
    const nucleo = b.dependency("nucleo", .{ .target = target, .optimize = optimize }).module("nucleo");
    const mod = b.addModule("kerf", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{ .{ .name = "bitdp", .module = bitdp }, .{ .name = "nucleo", .module = nucleo } },
    });
    const tests = b.addTest(.{ .root_module = mod });
    b.step("test", "Run tests").dependOn(&b.addRunArtifact(tests).step);

    const demo = b.addExecutable(.{
        .name = "kerf-demo",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/demo.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "kerf", .module = mod },
                .{ .name = "nucleo", .module = nucleo },
            },
        }),
    });
    b.step("demo", "Build the population off-target demo").dependOn(&b.addInstallArtifact(demo, .{}).step);
}

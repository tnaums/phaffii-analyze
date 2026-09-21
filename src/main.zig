const std = @import("std");
const fasta = @import("fasta");
const Io = std.Io;

const Komagataella = struct {
    plasmid: *fasta.DNA,
};

fn isInducible(k: Komagataella) bool {
    const d = std.mem.indexOf(u8, k.plasmid.sequence, "AGATCTAACATCCAAA");
    const e = std.mem.indexOf(u8, k.plasmid.sequence, "ATTCGAAACGA") orelse 0;

    if (d) |value| {
        if (value == 0 and e == 930) {
            return true;
        }
    }
    return false;
}

fn isSecreted(k: Komagataella) bool {
    if (std.mem.eql(u8, k.plasmid.sequence[940..949], "ATGAGATTT") and
            (std.mem.eql(u8, k.plasmid.sequence[1198..1207], "GCTGAAGCT"))
        ) {
        std.debug.print("{s}\n", .{k.plasmid.sequence[940..949]});
        std.debug.print("{s}\n", .{k.plasmid.sequence[1198..1207]});
        return true;
    }
    return false;
}

fn matureProtein(k: Komagataella, allocator: std.mem.Allocator) !void {
    const endpoint = std.mem.indexOf(u8, k.plasmid.sequence, "TCTAGA");
    if (endpoint) |v| {
        const sequence = k.plasmid.sequence[1207..v];
        var coding = try fasta.DNA.init(allocator, "mature", sequence);
        defer coding.deinit(allocator);
        try coding.addTranslation(allocator);
        var protein: fasta.Protein = undefined;
        defer protein.deinit(allocator);
        if (coding.translation) |value| {
            protein = try fasta.Protein.init(allocator, "mature", value[0]);
        }
        std.debug.print("{f}\n", .{protein});
        return;
    }
    std.debug.print("endpoint not found\n", .{});

}

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();
    
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.debug.print("Usage: {s} <filename.fa>\n", .{args[0]});
        return;
    }

    const filepath = args[1];
    const file = try std.Io.Dir.cwd().openFile(init.io, filepath, .{});
    defer file.close(init.io);

    // Starting DNA queue
    var queue: std.Io.Queue(fasta.DNA) = .init(&.{});
    var producer_task = try init.io.concurrent(fasta.parseDNA, .{ init.io, init.gpa, &queue, file });
    defer producer_task.cancel(init.io) catch {};

    var myPlasmid = try queue.getOne(init.io);
    try myPlasmid.addTranslation(init.gpa);
    defer myPlasmid.deinit(init.gpa);

    const k = Komagataella{
        .plasmid = &myPlasmid,
    };

    const fs = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k.plasmid});
    defer init.gpa.free(fs);
    try stdout.writeStreamingAll(init.io, fs);

    try myPlasmid.mapDNA(init.gpa, init.io, stdout);
    const b = isInducible(k);
    if (b) {
        try stdout.writeStreamingAll(init.io, "The plasmid is inducible!\n");
    }
    const b2 = isSecreted(k);
    if (b2) {
        try stdout.writeStreamingAll(init.io, "The protein is secreted.\n");
    }
    try matureProtein(k, init.gpa);
    
}

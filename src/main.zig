const std = @import("std");
const fasta = @import("fasta");
const Io = std.Io;

const promoter = enum { aox1, gap, unknown };
const secretion = enum { alpha, ost, cytoplasmic };

const Komagataella = struct {
    plasmid: fasta.DNA,
    promoter: promoter,
    coding: fasta.DNA,
    secretion: secretion,
    recombinant: fasta.Protein,

    pub fn init(io: Io, allocator: std.mem.Allocator, filepath: []const u8) !Komagataella {
        var dna = try fasta.parseSingleDNA(io, allocator, filepath);
        try dna.addTranslation(allocator);
        const p = checkPromoter(dna);
        var codingRegion = try getCodingRegion(allocator, dna, p);
        try codingRegion.addTranslation(allocator);
        const s = determineSecretion(codingRegion);
        const r = try recombinantProtein(allocator, codingRegion, s);
        return Komagataella{
            .plasmid = dna,
            .promoter = p,
            .coding = codingRegion,
            .secretion = s,
            .recombinant = r,
        };
    }

    pub fn deinit(self: Komagataella, allocator: std.mem.Allocator) void {
        self.plasmid.deinit(allocator);
        self.coding.deinit(allocator);
        self.recombinant.deinit(allocator);
    }

    fn checkPromoter(dna: fasta.DNA) promoter {
        const aox1Start = std.mem.indexOf(u8, dna.sequence, "AGATCTAACATCCAAA");
        const aox1End = std.mem.indexOf(u8, dna.sequence, "ATTCGAAACGA") orelse 0;

        const gapStart = std.mem.indexOf(u8, dna.sequence, "AGATCTTTTTTGTAGAAATG");
        const gapEnd = std.mem.indexOf(u8, dna.sequence, "TATTTCAATCAATTGAACAAC") orelse 0;

        if (aox1Start) |value| {
            if (value == 0 and aox1End == 930) {
                return promoter.aox1;
            }
        }
        if (gapStart) |value| {
            if (value == 0 and gapEnd == 459) {
                return promoter.gap;
            }
        }
        return promoter.unknown;
    }

    fn getCodingRegion(allocator: std.mem.Allocator, dna: fasta.DNA, p: promoter) !fasta.DNA {
        const newHeader = try std.fmt.allocPrint(allocator, "{s}|{s}", .{ dna.header, "coding"});
        defer allocator.free(newHeader);
        
        switch(p) {
            .aox1 => {
                const codingEnd = std.mem.indexOf(u8, dna.sequence, "GTTTGTAGCCTTAGA") orelse dna.sequence.len;
                return try fasta.DNA.init(allocator, newHeader, dna.sequence[940..codingEnd]);
            },
            .gap => {
                const codingEnd = std.mem.indexOf(u8, dna.sequence, "GTTTTAGCCTTAGAC") orelse dna.sequence.len;                
                return try fasta.DNA.init(allocator, newHeader, dna.sequence[492..codingEnd]);
            },
            .unknown => {unreachable;},
        }
    }

    fn determineSecretion(dna: fasta.DNA) secretion {
        if  (std.mem.startsWith(u8, dna.sequence, "ATGAGATTTCCTTCA")) {
            return secretion.alpha;
        }
        if (std.mem.startsWith(u8, dna.sequence, "ATGAGGCAGGTTTGG")) {
            return secretion.ost;
        }
        return secretion.cytoplasmic;
    }

    fn recombinantProtein(allocator: std.mem.Allocator, codingRegion: fasta.DNA, s: secretion) !fasta.Protein {
        const newHeader = try std.fmt.allocPrint(allocator, "{s}|{s}", .{ codingRegion.header[0..codingRegion.header.len - 7], "mature"});
        defer allocator.free(newHeader);
        var fullP: []u8 = undefined;

        if (codingRegion.translation) |value| {
            fullP = value[0];
        }
        var parts = std.mem.tokenizeScalar(u8, fullP, '*');
        var shortLeft = parts.next() orelse fullP;

        switch (s) {
            .alpha => {
                return try fasta.Protein.init(allocator, newHeader, shortLeft[89..]);
            },

            .ost => {
                return try fasta.Protein.init(allocator, newHeader, shortLeft[92..]);
            },
            .cytoplasmic => {
                return try fasta.Protein.init(allocator, newHeader, shortLeft);
            },
        }
    }

    pub fn format(self: Komagataella, writer: *Io.Writer) !void {
        const bp = self.plasmid.sequence.len;
        try writer.print("  plasmid: {s}\n   length: {d}bp\n", .{ self.plasmid.header, bp });
        switch (self.promoter) {
            .aox1 => {
                try writer.print(" promoter: methanol inducible aox1\n", .{});
            },
            .gap => {
                try writer.print(" promoter: constitutive gap promoter\n", .{});
            },
            .unknown => {
                try writer.print(" promoter: unknown, probably not a pPICZ expression plasmid\n", .{});
            },
        }
        switch (self.secretion) {
            .alpha => {
                try writer.print("secretion: α factor\n", .{});
            },
            .ost => {
                try writer.print("secretion: ost1\n", .{});
            },
            .cytoplasmic => {
                try writer.print("secretion: standard SSS not detected, may be cytoplasmic\n", .{});
            },
        }
        try writer.print("\n", .{});
        try writer.print("{f}\n", .{self.recombinant});
    }
};

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();

    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.debug.print("Usage: {s} <filename.fa>\n", .{args[0]});
        return;
    }

    const filepath = args[1];

    var k = try Komagataella.init(init.io, init.gpa, filepath);
    defer k.deinit(init.gpa);

    const fs = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k.coding});
    defer init.gpa.free(fs);
    try stdout.writeStreamingAll(init.io, fs);

    try k.coding.mapDNA(init.gpa, init.io, stdout, fasta.readingframe.first);

    try stdout.writeStreamingAll(init.io, "----------------------------------------\n");

    const kf = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k});
    defer init.gpa.free(kf);
    try stdout.writeStreamingAll(init.io, kf);

}
